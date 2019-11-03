#! /usr/bin/env python

"""
Michael Tucker 2018
Simple script for extracting and combining SNIFS spectra
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from astropy.time import Time
import sys
import argparse
import ipdb
from fluxcalSNIFS import fluxcalSNIFS

try:
	from spectres import spectres
except ImportError:
	print('Could not import SpectRes! Need to install?')
	sys.exit(1)


BAD_REGIONS = [[4900.0,5200.0], [7590.0, 7710.0]]

def main(directory='./', binsz=7.0, plot=True, overwrite=False, blue=True, red=True, fluxcal=False, 
									std_coord_match=False, verbose=True, names=None, avoid=None, combine='mask'):
	if not os.path.exists(directory): raise IOError('Could not find directory %s' % directory)
	if not os.path.isdir(directory): raise IOError('%s is not a directory' % directory)
	if binsz <= 0.0: raise ValueError('Provided bin size is <= 0!')


	if verbose:
		print('-'*30)
		print('Running SNIFS extraction script')
		print('-'*30)

	fnames = sorted(FindFiles(directory, verbose))


	written_files = []
	RSPEX = []
	BSPEX = []
	
	for fbase in fnames:
#		ipdb.set_trace()
		if verbose: print('\n\n'+'-'*30+'\n'+'Processing base %s...' % fbase+'\n'+'-'*30)
		if red:
			rff = os.path.join(directory,fbase + 'R.fits')
			if not os.path.exists(rff):
				print('WARNING: Could not find R spectrum (%s), skipping' %rff)
				rspec = None
				rff = None
			else:
				rspec = Spectrum(rff)
				if len(avoid) > 0:
					if rspec.obj in avoid: continue
				if len(names) > 0:
					if rspec.obj not in names: continue

				rvar = os.path.join(directory,'var_'+fbase+'R.fits')
				if not os.path.exists(rvar): raise IOError('Could not find variance file for %s' % rff)

				rspec.add_err(rvar)
				RSPEX.append(rspec)
		else:
			rspec = None

		if blue:
			bff = os.path.join(directory,fbase + 'B.fits')
			if not os.path.exists(bff):
				print('WARNING: Could not find B spectrum (%s), skipping' %bff)
				bspec = None
				bff = None
			if bff != None:
				bspec = Spectrum(bff)
				if len(avoid) > 0:
					if bspec.obj in avoid: continue
				if len(names) > 0:
					if bspec.obj not in names: continue

				bvar = os.path.join(directory,'var_'+fbase+'B.fits')
				if not os.path.exists(bvar): raise IOError('Could not find variance file %s' % bvar)

				bspec.add_err(bvar)
				BSPEX.append(bspec)
		else:
			bspec = None

	if fluxcal:
		if red: fluxcalSNIFS(RSPEX, plot, std_coord_match, verbose)
		if blue: fluxcalSNIFS(BSPEX, plot, std_coord_match, verbose)

	if red and blue:
		bnames = [spec.object for spec in BSPEX]
		rnames = [spec.object for spec in RSPEX]
		done_names = []
		spex = []
		for name,spec in zip(bnames, BSPEX):
			if not name in rnames: spex.append(spec)
			else: 
				jj = rnames.index(name)
				fullspec = CombineSpecArms(spec, RSPEX[jj], binsz=binsz, combine=combine)
				spex.append(fullspec)
			done_names.append(name)
		if not all([True if name in done_names else False for name in rnames]):
			for name,spec in zip(rnames, RSPEX):
				if name in done_names: continue
				spex.append(spec)
				done_names.append(name)
	elif red:
		spex = RSPEX.copy()
	elif blue:
		spex = BSPEX.copy()

	written_files = []
	for spec in spex:
		WriteSpectrum(spec, written_files, plot=plot, verbose=verbose, overwrite=overwrite)

def CombineSpecArms(bspec=None, rspec=None, combine='mask', binsz=7.0, verbose=True):
		if rspec and bspec:
			if verbose: print('\tCombing blue and red arms.')
			spec = CombineSpectra(bspec, rspec, combine=combine, rebin=binsz)
		elif bspec:
			if verbose: print('\tWARNING: Only blue arm found.')
			spec = bspec
		elif rspec:
			if verbose: print('\tWARNING: Only red arm found.')
			spec = rspec
		else:
			raise RuntimeError('Something went wrong processing spectra...')
		return spec


def WriteSpectrum(spec, written_files, plot=True, verbose=True, overwrite=False):
	if spec == None: return
	newff = spec.object +'-'+str(round(spec.mjd, 4))+'-SNIFS-0.dat'

	count = 0
	while newff in written_files:
		newff = np.str.replace(newff, str(count)+'.dat', str(count+1)+'.dat')
		count += 1

	if os.path.exists(newff) and not overwrite: 
		print('%s already exists! See overwrite flags, skipping...' % newff)
		return

	if plot:
		plt.plot(spec.wl, spec.fl, 'k-', drawstyle='steps-mid')
		plt.title(spec.object)
		plt.show()

	if spec.fluxcal:
		FCAL = 'ABSOLUTE'
	else:
		FCAL = 'RELATIVE'

	if verbose: print('\tWriting spectrum to %s' % newff)
	with open(newff, 'w') as ofile:
		ofile.write('#OBJECT = %s\n' % spec.object)
		ofile.write('#MJD-OBS = %lf\n' % spec.mjd)
		ofile.write('#DATE-OBS = %s\n' % spec.date)
		ofile.write('#EXPTIME = %lf\n' % spec.expt)
		ofile.write('#AIRMASS = %lf\n# \n' % spec.am)
		ofile.write('#TELESCOPE = UH88\n')
		ofile.write('#INSTRUMENT = SNIFS\n')
		ofile.write('#FLUXCAL = %s\n' % FCAL)
#		ofile.write('#BLUE_FNAME = %s\n' % spec.bfluxf)
#		ofile.write('#BLUE_VARF = %s\n' % spec.bvarf)
#		ofile.write('#RED_FNAME = %s\n' % spec.rfluxf)
#		ofile.write('#RED_VARF = %s\n# \n' % spec.rvarf)		
		ofile.write('#WL FL ERR\n')
		ofile.write('#[Ang] [erg/s/cm2/A] [erg/s/cm2/A]\n')
		for item in zip(spec.wl, spec.fl, spec.err): ofile.write('%lf %le %le\n' % item)

	written_files.append(newff)



def FindFiles(directory, verbose):
	fnames = list(set([ff[:-6] for ff in os.listdir(directory) if ff.startswith('spec_') and ff.endswith('.fits')]))
	if verbose: print('\t Found %d SNIFS spectra.' % len(fnames))
	if len(fnames) == 0: print('No spec*.fits files found in dir: %s,\n\t maybe need -d flag?' % directory)
	return fnames



class Spectrum():
	def __init__(self, fitsfile):
		self.fname = fitsfile
		self.fl, header = fits.getdata(fitsfile, header=True)
		self.object = header['OBJECT']

		self.date=header['DATE-OBS']

		if '/' in self.date:
			self.date = np.str.replace(date, '/', '-')

		self.mjd = Time(self.date).mjd
		self.expt = header['EXPTIME']
		self.am = header['AIRMASS']
		self.instrument='SNIFS'
		self.telescope='UH88'

		rastr = header['OBJRA'].split(':')
		self.ra = float(rastr[0])*15.0 + float(rastr[1])/60.0 + float(rastr[2])/3600.0
		if header['OBJDEC'][0] == '-':
			sign = -1.0
		else:
			sign = 1.0
		decstr = header['OBJDEC'].split(':')
		self.dec = sign * ( float(decstr[0]) + float(decstr[1])/60.0 + float(decstr[2])/3600.0 )

		self.wl = header['CRVAL1'] + header['CDELT1']*np.arange(len(self.fl))
		self.fluxcal = False

	def add_err(self, fitsfile, verify=True):
		self.err_file = fitsfile
		err, header = fits.getdata(fitsfile, header=True)

		if verify:
			date = np.str.replace(header['DATE-OBS'], '/', '-')
			mjd = Time(date).mjd
			assert (mjd - self.mjd)*24.0*60.0 < 1.0 
			assert header['EXPTIME'] == self.expt and header['AIRMASS'] == self.am and header['OBJECT'] == self.object
		self.err = np.sqrt(err)

	def GetExposure(self):
		return int(self.fname.split('_')[3])

	def rebin(self, binsz):
		newwl = np.arange(self.wl.min()+2.0*binsz, self.wl.max()-2.0*binsz, binsz)
		newfl, newerr = spectres(newwl, self.wl, self.fl, self.err)
		self.wl = newwl
		self.fl = newfl
		self.err = newerr

def CombineSpectra(bspec, rspec, combine='wavg', rebin=2.0):

	if combine =='smart':
		full_wave, full_fl, ferr = smartCombine(bspec=bspec, rspec=rspec, rebin=rebin)

	elif combine in ['med', 'avg', 'wavg']:
		wl0 = bspec.wl.min()
		wl1 = rspec.wl.max()
		full_wave = np.arange(wl0+rebin, wl1-rebin, rebin)
		bfl = np.interp(full_wave, bspec.wl, bspec.fl, left=0.0, right=0.0)
		berr = np.interp(full_wave, bspec.wl, bspec.err, left=0.0, right=0.0)
		rfl = np.interp(full_wave, rspec.wl, rspec.fl, left=0.0, right=0.0)
		rerr = np.interp(full_wave, rspec.wl, rspec.err, left=0.0, right=0.0)

		if combine == 'med':
			full_fl = np.median(np.stack([bfl, rfl]), axis=0)
		elif combine == 'avg':
			full_fl = np.average(np.stack([bfl, rfl]), axis=0)
		elif combine=='wavg':
			full_fl = np.average(np.stack([bfl, rfl]), axis=0, weights=np.stack([berr**-2.0, rerr**-2.0]))
		else:
			raise RuntimeError('Something went wrong in CombineSpectra...')

		ipdb.set_trace()
		ferr = np.sqrt(berr**2.0 + rerr**2.0)

		

	elif combine=='mask':
		bkeep = np.where(bspec.wl <= BAD_REGIONS[0][0])[0]
		rkeep = np.where(rspec.wl >= BAD_REGIONS[0][1])[0]
		wl = np.append(bspec.wl[bkeep], rspec.wl[rkeep])
		fl = np.append(bspec.fl[bkeep], rspec.fl[rkeep])
		fe = np.append(bspec.err[bkeep], rspec.err[rkeep])

		wl0 = bspec.wl.min()
		wl1 = rspec.wl.max()
		full_wave = np.arange(wl0+rebin, wl1-rebin, rebin)
		full_fl, ferr = spectres(full_wave, wl, fl, fe)

	elif combine =='snid':
		bkeep = np.where(bspec.wl <= BAD_REGIONS[0][0])[0]
		rkeep = np.where(rspec.wl <= BAD_REGIONS[0][1])[0]
		full_wave = np.arange(bspec.wl.min()+binsz, rspec.wl.max()-binsz, binsz)
		full_fl = np.interp(full_wave, )

	else:
		raise ValueError('Unknown combine type %s' % combine)


	spec = bspec
	spec.wl = full_wave
	spec.fl = full_fl
	spec.err = ferr

	spec.bfluxf = bspec.fname
	spec.bvarf = bspec.err_file
	spec.rfluxf = rspec.fname
	spec.rvarf = rspec.err_file


	return spec

def smartCombine(bspec, rspec, rebin=5.0):
	minwl = bspec.wl.min()
	bedge = bspec.wl.max()
	maxwl = rspec.wl.max()
	redge = rspec.wl.min()

	newwl = np.arange(minwl+rebin, maxwl-rebin, rebin)
	newfl = np.zeros_like(newwl)
	newerr = np.zeros_like(newwl)

	indices = np.arange(newwl.shape[0], dtype=int)
	blue_only = indices[newwl <= redge+0.5*rebin]
	newfl[blue_only], newerr[blue_only] = spectres(newwl[blue_only], bspec.wl, bspec.fl, bspec.err)

	red_only = indices[newwl >= bedge-0.5*rebin]
	newfl[red_only], newerr[red_only] = spectres(newwl[red_only], rspec.wl, rspec.fl, rspec.err)

	overlap = np.where(newfl == 0.)[0]
	bfl_olap, berr_olap = spectres(newwl[overlap], bspec.wl, bspec.fl, bspec.err)
	rfl_olap, rerr_olap = spectres(newwl[overlap], rspec.wl, rspec.fl, rspec.err)

	bweights = berr_olap**-2. * np.linspace(1, 0, len(overlap))
	rweights = rerr_olap**-2. * np.linspace(0, 1, len(overlap))

	newfl_olap = []
	newerr_olap = []
	for bf, be, bw, rf, re, rw in zip(bfl_olap, berr_olap, bweights, rfl_olap, rerr_olap, rweights):
		avg, std = weighted_avg_and_std(np.array([bf, rf]), np.array([bw, rw]))
		newfl_olap.append(avg)
		newerr_olap.append(std)

	ipdb.set_trace()

	newfl[overlap] = np.array(newfl_olap)
	newerr[overlap] = np.sqrt(berr_olap**2.0 + rerr_olap**2. + (0.1*newfl[overlap])**2. )
	return newwl, newfl, newerr


def weighted_avg_and_std(values, weights):
	average = np.average(values, weights=weights)
	# Fast and numerically precise:
	variance = np.average((values-average)**2, weights=weights)
	return (average, np.sqrt(variance))


if __name__=='__main__':
	description="Script for extracting SNIFS spectra already reduced by the pipeline."
	usage = "./extractSNIFS.py [options]"
	parser = argparse.ArgumentParser(description=description, usage=usage)
	parser.add_argument('-d', '--dir', help='Directory where spec_*.fits files are located. Default: CWD', default=os.getcwd(), type=str)
	parser.add_argument('-w', '--binsz', help='Output bin size of spectrum: default: 7.0 (approx spec resolution)', default=7.0, type=float)
	parser.add_argument('-p', '--plot', help='Turn on (1)/off(0) plotting. Default: 1 (on)', default=1, type=int, choices=[0,1])	
	parser.add_argument('-b', '--blue', help='Extract blue channel? (0/1) Default: 1', default=1, type=int, choices=[0,1])
	parser.add_argument('-r', '--red', help='Extract red channel? (0/1) Default: 1', default=1, type=int, choices=[0,1])
	parser.add_argument('-q', '--quiet', help='supress output? Default: False', default=False, action='store_true')
	parser.add_argument('-n', '--name', help='Object names to be extracted, can be used multiple times: -n NAME1 -n NAME2 -- or -- -n NAME1 NAME2 etc. Default: None (extract all)',
		default='', type=str)
	parser.add_argument('-a', '--avoid', help='Objects to avoid (not extract). Same as above, can be re-used for multiple objects.', default='', type=str)
	parser.add_argument('--coord_match', help='Coordinate-match stanbdard stars instead of name-match. Default: False', default=False, action='store_true')
	parser.add_argument('-o', '--overwrite', help='Overwrite existing extraced spectra? Default: False', default=False, action='store_true')
	parser.add_argument('-c', '--combine', help='Combine method for blue/red arms. Default: mask', default='mask', 
		choices=['med', 'avg', 'wavg', 'mask', 'snid', 'smart'], type=str)
	parser.add_argument('--fluxcal', help='Apply flux calibrations from standard stars? Default: False', action='store_true', default=False)
	args = parser.parse_args()
	main(args.dir, args.binsz, args.plot, args.overwrite, args.blue, args.red, args.fluxcal, args.coord_match, not args.quiet, args.name, args.avoid,args.combine)
