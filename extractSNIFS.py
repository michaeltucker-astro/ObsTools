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
from spectra import CombineSpec

try:
	from SpectRes import spectres
except ImportError:
	print('Could not import SpectRes! Need to install?')
	sys.exit(1)


BAD_REGIONS = [[4900.0,5200.0], [7590.0, 7710.0]]

def main(directory, binsz, plot, overwrite, mask, blue, red, verbose, names, avoid, combine):
	if not os.path.exists(directory): raise IOError('Could not find directory %s' % directory)
	if not os.path.isdir(directory): raise IOError('%s is not a directory' % directory)
	if binsz <= 0.0: raise ValueError('Provided bin size is <= 0!')


	if verbose:
		print('-'*30)
		print('Running SNIFS extraction script')
		print('-'*30)

	fnames = FindFiles(directory, verbose)
	for fbase in fnames:
#		ipdb.set_trace()
		if verbose: print('\n\n'+'-'*30+'\n'+'Processing base %s...' % fbase+'\n'+'-'*30)
		if red:
			rff = os.path.join(directory,fbase + 'R.fits')
			if not os.path.exists(rff):
				print('WARNING: Could not find R spectrum (%s), skipping' %rff)
				rspec = None
				rff = None
			if rff != None:
				rspec = Spectrum(rff)
				if len(avoid) > 0:
					if rspec.obj in avoid: continue
				if len(names) > 0:
					if rspec.obj not in names: continue

				rvar = os.path.join(directory,'var_'+fbase+'R.fits')
				if not os.path.exists(rvar): raise IOError('Could not find variance file for %s' % rff)

				rspec.add_err(rvar)
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
		else:
			bspec = None


		if rspec and bspec:
			if verbose: print('\tCombing blue and red arms.')
			spec = CombineSpectra(bspec, rspec, combine=combine)
		elif bspec:
			if verbose: print('\tWARNING: Only blue arm found.')
			spec = bspec
		elif rspec:
			if verbose: print('\tWARNING: Only red arm found.')
			spec = rspec
		else:
			raise RuntimeError('Something went wrong processing spectra...')

		spec.rebin(binsz)

		newff = spec.obj +'-'+str(round(spec.mjd, 4))+'-SNIFS.dat'
		if os.path.exists(newff) and not overwrite: 
			print('%s already exists! See overwrite flags, skipping...' % newff)
			continue


		if plot:
			plt.plot(spec.wl, spec.fl, 'k-', drawstyle='steps-mid')
			plt.title(spec.obj)
			plt.show()

		if verbose: print('\tWriting spectrum to %s' % newff)
		with open(newff, 'w') as ofile:
			ofile.write('#OBJECT = %s\n' % spec.obj)
			ofile.write('#MJD-OBS = %lf\n' % spec.mjd)
			ofile.write('#DATE-OBS = %s\n' % spec.date)
			ofile.write('#EXPTIME = %lf\n' % spec.expt)
			ofile.write('#AIRMASS = %lf\n' % spec.am)
			ofile.write('#TELESCOPE = UH88\n')
			ofile.write('#INSTRUMENT = SNIFS\n')
			ofile.write('#WL FL ERR\n')
			ofile.write('#[Ang] [erg/s/cm2/A] [erg/s/cm2/A]\n')
			for item in zip(spec.wl, spec.fl, spec.err): ofile.write('%lf %le %le\n' % item)



def FindFiles(directory, verbose):
	fnames = list(set([ff[:-6] for ff in os.listdir(directory) if ff.startswith('spec_') and ff.endswith('.fits')]))
	if verbose: print('\t Found %d SNIFS spectra.' % len(fnames))
	return fnames



def ExtractSNIFS(directory=None, binsz=1.0,plot=True,overwrite=False,verbose=True, mask_dichroic=False, mask_Aband=False, blue=True, red=True):
	if directory == None: directory = os.getcwd()
	fnames = [os.path.join(directory,ff) for ff in os.listdir(directory) if ff.startswith('spec') and ff.endswith('.fits')]
	if verbose: print('Found %d files for extraction' % len(fnames))

	varnames = []
	for ff in (fnames):
		if ff.endswith('R.fits'): 
			rff = ff
			bff = np.str.replace(rff, 'R', 'B')
		elif ff.endswith('B.fits'):
			bff = ff
			rff = np.str.replace(bff, 'B', 'R')

		varbff = os.path.join(directory, 'var_'+os.path.split(bff)[-1])
		varrff = os.path.join(directory, 'var_'+os.path.split(rff)[-1])
		
		if blue:
			if not os.path.exists(bff): 
				if skipbad: continue
				else: raise IOError('Could not find file for %s' % bff)			
			if not os.path.exists(varbff): 
				if skipbad: continue
				else: raise IOError('Could not find file for %s' % varbff)
	
			#bwl, bfl, berr, bobj, bdate, bexpt, bam = read_spec(bff, varbff)
			bspec = Spectrum(bff)
			bspec.add_err(varbff)
		if red:			
			if not os.path.exists(rff): 
				if skipbad: continue
				else: raise IOError('Could not find file for %s' % rff)
			if not os.path.exists(varrff): 
				if skipbad: continue
				else: raise IOError('Could not find file for %s' % varrff)

			#rwl, rfl, rerr, robj, rdate, rexpt, ram = read_spec(rff, varrff)
			rspec = Spectrum(rff)
			rspec.add_err(varrff)
		if blue and red:
			obj = bspec.obj
			expt = np.mean([bspec.expt, rspec.expt])
			am = np.mean([bspec.am, rspec.am])
			date = bspec.date

			wl, fl = CombineSpec([bspec.wl, rspec.wl], [bspec.fl, rspec.fl], rebin=binsz)
			wl, err = CombineSpec([bspec.wl, rspec.wl], [bspec.err, rspec.err], rebin=binsz)
			spec = bspec
			spec.wl = wl
			spec.fl =fl
			spec.err = err

		elif blue:
			obj = bspec.obj
			bwl1, bfl = CombineSpec([bspec.wl, bspec.wl], [bspec.fl, bspec.fl], rebin=binsz)
			bwl1, berr = CombineSpec([bspec.wl, bspec.wl], [bspec.err, bspec.err], rebin=binsz)			
			spec = bspec
			spec.wl = bwl1
			spec.fl = bfl
			spec.err = berr


		elif red:
			obj = rspec.obj
			rwl1, rfl = CombineSpec([rspec.wl, rspec.wl], [rspec.fl, rspec.fl], rebin=binsz)
			rwl1, rerr = CombineSpec([rspec.wl, rspec.wl], [rspec.fl, rspec.fl], rebin=binsz)
			spec = rspec
			spec.wl = rwl1
			spec.fl = rfl
			spec.err = rerr
		else:
			raise RuntimeError('Something went wrong in ExtractSNIFS!')

		newpath = obj+'-'+str(round(spec.mjd,4))+'-SNIFS.dat'
		if os.path.exists(newpath) and not overwrite: continue
		elif os.path.exists(newpath) and overwrite:
			resp = input("Filename %s already exists, overwrite? [y/N] >" % newpath)
			if resp.lower().strip() != 'y': continue

		if verbose: print ('\nExtracting spectra for %s' % obj)

		if plot:
			plt.plot(spec.wl, spec.fl, color='b', drawstyle='steps-mid', ls='-', lw=1.0)
			for reg in BAD_REGIONS:
				plt.fill_betweenx([np.min(spec.fl),np.max(spec.fl)], reg[0], reg[1], facecolor='k', alpha=0.5)
			plt.title(obj)
			plt.xlabel('Wavelength (A)')
			plt.ylabel('Flux')
			plt.xlim(spec.wl.min(), spec.wl.max())
			plt.show()

		if mask_dichroic:
			keep = np.where((spec.wl <= BAD_REGIONS[0][0])|(spec.wl >= BAD_REGIONS[0][1]))[0]
			spec.wl = spec.wl[keep]
			spec.fl = spec.fl[keep]
			spec.err = spec.err[keep]

		if mask_Aband and red:
			keep = np.where((spec.wl <= BAD_REGIONS[1][0])|(spec.wl >= BAD_REGIONS[1][1]))[0]
			spec.wl = spec.wl[keep]
			spec.fl = spec.fl[keep]
			spec.err = spec.err[keep]

		with open(newpath, 'w') as ofile:
			ofile.write('# MJD-OBS = %lf\n' % spec.mjd)
			ofile.write('# DATE-OBS = %s\n' % spec.date)
			ofile.write('# EXPTIME = %lf\n' % spec.expt)
			ofile.write('# AIRMASS = %4.3f\n' % spec.am)
			ofile.write('# TELESCOPE = UH88\n# INSTRUMENT = SNIFS\n')
			for w, f, e in zip(spec.wl, spec.fl, spec.err): ofile.write('%le %le %le\n' % (w,f, e))
		print ('\tSpectrum written to ./%s' % newpath)

		for ff in [bff, rff, varbff, varrff]: 
			try:
				fnames.remove(ff)
			except:
				continue

class Spectrum():
	def __init__(self, fitsfile):
		self.fname = fitsfile
		self.fl, header = fits.getdata(fitsfile, header=True)
		self.obj = header['OBJECT']

		self.date=header['DATE-OBS']

		if '/' in self.date:
			self.date = np.str.replace(date, '/', '-')

		self.mjd = Time(self.date).mjd
		self.expt = header['EXPTIME']
		self.am = header['AIRMASS']
		self.instrument='SNIFS'
		self.telescope='UH88'

		self.wl = header['CRVAL1'] + header['CDELT1']*np.arange(len(self.fl))

	def add_err(self, fitsfile, verify=True):
		err, header = fits.getdata(fitsfile, header=True)
		if verify:
			date = np.str.replace(header['DATE-OBS'], '/', '-')
			mjd = Time(date).mjd
			assert (mjd - self.mjd)*24.0*60.0 < 1.0 
			assert header['EXPTIME'] == self.expt and header['AIRMASS'] == self.am and header['OBJECT'] == self.obj
		self.err = err

	def GetExposure(self):
		return int(self.fname.split('_')[3])

	def rebin(self, binsz):
		newwl = np.arange(self.wl.min()+2.0*binsz, self.wl.max()-2.0*binsz, binsz)
		newfl, newerr = spectres.spectres(self.wl, self.fl, newwl, self.err)
		self.wl = newwl
		self.fl = newfl
		self.err = newerr

def CombineSpectra(bspec, rspec, combine='wavg', rebin=2.0):

	if combine in ['med', 'avg', 'wavg']:
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
		full_wave = np.append(bspec.wl[bkeep], rspec.wl[rkeep])
		full_fl = np.append(bspec.fl[bkeep], rspec.fl[rkeep])
		ferr = np.append(bspec.err[bkeep], rspec.err[rkeep])
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

	return spec

if __name__=='__main__':
	description="Script for extracting SNIFS spectra already reduced by the pipeline."
	parser = argparse.ArgumentParser(description)
	parser.add_argument('-d', '--dir', help='Directory where spec_*.fits files are located. Default: CWD', default=os.getcwd(), type=str)
	parser.add_argument('-w', '--binsz', help='Output bin size of spectrum: default: 2.0', default=2.0, type=float)
	parser.add_argument('-p', '--plot', help='Turn on (1)/off(0) plotting. Default: 1 (on)', default=1, type=int, choices=[0,1])
	parser.add_argument('-m', '--mask', help='Whether to mask the dichroic, the A-band, both, or neither. Default: both', default='both', type=str, 
		choices=['both', 'neither', 'dichroic', 'Aband'])
	parser.add_argument('-b', '--blue', help='Extract blue channel? (0/1) Default: 1', default=1, type=int, choices=[0,1])
	parser.add_argument('-r', '--red', help='Extract red channel? (0/1) Default: 1', default=1, type=int, choices=[0,1])
	parser.add_argument('-v', '--verbose', help='Verbose output off/on (0/1) Default: 0', default=0, type=int, choices=[0,1])
	parser.add_argument('-n', '--name', help='Object names to be extracted, can be used multiple times: -n NAME1 -n NAME2 -- or -- -n NAME1 NAME2 etc. Default: None (extract all)',
		default='', type=str)
	parser.add_argument('-a', '--avoid', help='Objects to avoid (not extract). Same as above, can be re-used for multiple objects.', default='', type=str)
	parser.add_argument('-o', '--overwrite', help='Overwrite existing extraced spectra? [0,1] Default: 0', default=0, choices=[0,1], type=int)
	parser.add_argument('-c', '--combine', help='Combine method for blue/red arms. Default: waverage', default='wavg', 
		choices=['med', 'avg', 'wavg', 'mask', 'snid'], type=str)

	args = parser.parse_args()
	main(args.dir, args.binsz, args.plot, args.overwrite, args.mask, args.blue, args.red, args.verbose, args.name, args.avoid,args.combine)
