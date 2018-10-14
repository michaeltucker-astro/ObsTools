"""
Michael Tucker 2018
Simple tools to make observing easier

"""
import argparse
from datetime import date, timedelta
import requests
import sys,os
from pandas import read_csv, DataFrame
from io import StringIO
from astropy.coordinates import SkyCoord
import numpy as np
from collections import OrderedDict


def CSVtoJsky(fname, outdir=None):
	"""
	function: CSVtoJsky

	Converts a CSV file to a file parse-able by JSkyCalc
	CSV file must have colnames Target, RA, DEC

	Inputs:
		fname: filename of csv file
		outdir (optional): where to write the output file. Defaults to cwd

	Outputs:
		None
	"""

	def ConvertCoords(ra,dec):
		try:
			ra = float(ra) * 24.0 / 360.0
			hour = int(ra)
			ra = (ra - hour)*60.0
			minutes = int(ra)
			sec = round((ra - minutes) * 60.0, 2)
			rastr = '%02d %02d %05.2f' % (hour, minutes, sec)

			dec = float(dec)
			if dec >= 0.0: sign = '+'
			else: sign = '-'
			deg = int(dec)
			dec = (dec - deg) * 60.0
			dmin = int(dec)
			dsec = round((dec - dmin) * 60.0, 2)
			decstr = '%s%02d %02d %05.2f' % (sign,deg, dmin, dsec)

		except:
			if ':' in ra: ra = ra.split(':')
			else: ra = ra.split()
			rastr = '%02d %02d %05.2f' % (int(ra[0]), int(ra[1]), float(ra[2]))

			if ':' in dec: dec = dec.split(':')
			else: dec = dec.split()
			if int(dec[0]) >= 0.0: sign = '+'
			else: sign = '-'
			decstr = '%s%02d %02d %05.2f' % (sign,abs(int(dec[0])), int(dec[1]), float(dec[2]))
		finally:
			return rastr, decstr


	if not os.path.exists(fname):
		raise IOError('Cant find file %s!' % fname)

	if outdir == None: outdir = os.getcwd()+'/'


	data = read_csv(fname)
	if not all([True if colname in data.columns else False for colname in ['Object', 'RA', 'DEC']]):
		raise KeyError('CSV file must have colnames Object, RA and DEC!')


	ofile = open(outdir+'jsky.list','w')
	ofile.write('#List of targets for jskycalc\n')
	ofile.write('#File generated by CSVtoJsky\n')
	ofile.write('# Target \t RA \t\t DEC \t epoch\n#\n')

	line = '%s %s   %s  2000.0\n'
	for obj, ra, dec in zip(data['Object'], data['RA'], data['DEC']):
		if obj.startswith('#'): continue
		rastr, decstr = ConvertCoords(ra,dec)
		while len(obj) < 13: obj = obj + ' '
		ofile.write(line % (obj, rastr, decstr))

	ofile.close()
	print ('File written to %s' % (outdir+'jsky.list'))
	return

def FinderChart(ra, dec):
	import subprocess

	if ':' in ra:
		ra = ra.strip().split(':')
		ra = float(ra[0])*360.0/24.0 + float(ra[1])/60.0 + float(ra[2])/3600.0
	else:
		ra = float(ra.strip())

	if ':' in dec:
		if dec.startswith('-'): 
			dec = dec[1:]
			sign = -1.0
		else: sign = 1.0

		dec = dec.strip().split(':')
		dec = float(dec[0]) + float(dec[1])/60.0 + float(dec[2])/3600.0
		dec *= sign
	else:
		dec = float(dec.strip())

	baseURL = 'http://astro.subhashbose.com/render_AlL.php?RA=%lf\&DEC=%lf' % (ra,dec)
	subprocess.run('google-chrome --new-window %s' % baseURL, shell=True)

def Scheduler(fname, penalty=2.0):
	data = read_csv(fname)
	if not all([True if colname in data.columns else False for colname in ['Object','RA', 'DEC']]):
		raise KeyError('CSV file must have colnames Object, RA and DEC!')

	names = data['Object'].as_matrix()
	RA = data['RA'].as_matrix()
	DEC = data['DEC'].as_matrix()


	print('IMPLEMENTATION INCOMPLETE! exiting...')


def MakeTNSlist(ofile, declim, maglow, maghigh, Ndays, verbose):


	end = date.today()
	start = end - timedelta(days=Ndays)

	if verbose:
		print('#'*30)
		print('Querying TNS unclassified list')
		print('Parameters:')
		print('\tStart date: %s' % str(start))
		print('\tEnd date: %s' % str(end))
		print('\tNdays = %2d' % Ndays)
		print('\tMin mag: %3.1f' % maglow)
		print('\tMax mag: %3.1f' % maghigh)
		print('\tDec limit: %3.1f deg\n\n' % declim)

	baseURL = 'https://wis-tns.weizmann.ac.il/search?isTNS_AT=yes&public=all&unclassified_at=1&date_start%%5Bdate%%5D=%s&date_end%%5Bdate%%5D=%s&discovery_mag_min=%3.1f&discovery_mag_max=%3.1f' % (start, end, maglow, maghigh)
	endURL = '&num_page=1000&display%5Bredshift%5D=1&display%5Bhost_redshift%5D=1&display%5Binternal_name%5D=1&display%5Bdiscoverymag%5D=1&display%5Bdiscmagfilter%5D=1&display%5Bdiscoverydate%5D=1&format=csv' 
	fullURL = baseURL+endURL


	if verbose:print('Executing TNS request...')
	resp = requests.get(baseURL+endURL)
	if resp.status_code != 200:
		raise ValueError('Error retrieving list from TNS! Check internet?')
	else:
		if verbose: print('\tComplete!')

	if verbose: print('Parsing data...')
	table = read_csv(StringIO(resp.content.decode()))
	ATname = table['Name'].as_matrix()

	RA = [ra[0]+'h'+ra[1]+'m'+ra[2]+'s' for ra in [ra.split(':') for ra in table['RA'].as_matrix()]]
	DEC = [dec[0]+'d'+dec[1]+'m'+dec[2]+'s' for dec in [dec.split(':') for dec in table['DEC'].as_matrix()]]
	coords = SkyCoord(RA,DEC)
	RA = coords.ra.deg
	DEC = coords.dec.deg
	del coords

	redshift = table['Redshift'].as_matrix()
	hostz = table['Host Redshift'].as_matrix()
	name = table['Disc. Internal Name'].as_matrix()
	discmag = table['Discovery Mag'].as_matrix()
	discfilter = table['Discovery Mag Filter'].as_matrix()
	discdate = table['Discovery Date (UT)'].as_matrix()

	del table

	assert len(redshift[redshift > 0.0])==0
	del redshift

	keep = np.where(
		(DEC >= declim) &
		(discmag >= maglow) &
		(discmag <= maghigh)
		)[0]

	if verbose: print('\tComplete!')

	if verbose: print('Writing DataFrame...')

	df = DataFrame(OrderedDict({
		'Object':name[keep],
		'ATname':ATname[keep],
		'RA':RA[keep],
		'DEC':DEC[keep],
		'redshift':hostz[keep],
		'DiscMag':discmag[keep],
		'DiscFilter':discfilter[keep],
		'DiscDate':discdate[keep]
		}))


	if os.path.exists(ofile):
		print('-'*10)
		resp = 'null'
		while resp.strip().lower() not in ['','n','y'] and os.path.exists(ofile):
			resp = input('%s already exists, overwrite? [Y/n] >' % ofile).strip().lower()
			if resp in ['','y']: break
			elif resp == 'n':
				ofile = input('New output file (csv) >')

			else:
				print('Unrecognized input %s, ignoring...' % resp)

	df.to_csv(ofile, index=False)
	if verbose:print('DataFrame written to %s'%ofile)
	if verbose: print('total objects: %d' % len(keep))



if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Contains observing tools: CSVtoJsky, FinderChart')
	parser.add_argument('function', help='Which function to run', choices=['finderchart','fchart', 'fc','csv2jsky','c2j', 'scheduler', 'sched','tnslist', 'tns'], type=str)
	parser.add_argument('-r','--ra', help='RA for finderchart (deg)', type=str)
	parser.add_argument('-d','--dec', help='DEC for finderchart (deg)', type=str)
	parser.add_argument('-i', '--internet', help='Internet browser to use for finderchart. default: google-chrome', default='google-chrome', 
						choices=['google-chrome', 'firefox', 'safari'],type=str)

	parser.add_argument('-f','--fname', help='CSV file for CSVtoJsky or Scheduler', type=str)
	parser.add_argument('-p','--DECpenalty', help='Penalty factor for declination movements in Scheduler. Default: 2', default=2, type=float)
	parser.add_argument('-dl', '--dec-limit', help='Dec. limit for TNS target list. Default: -30 [deg]', default=-30.0, type=float)
	parser.add_argument('-mu', '--mag-upper', help='Mag upper limit for TNS target list. Default: 21 [mag]', default=21.0, type=float)
	parser.add_argument('-ml', '--mag-lower', help='Mag lower limit for TNS target list. Default: 10 [mag]', default=10.0, type=float)
	parser.add_argument('-o', '--output', help='Output filename for TNS target list. Default: tns-list.csv', default='tns-list.csv', type=str)
	parser.add_argument('-t', '--time', help='# of days to go back in TNS query. default: 7', default=7.0, type=float)
	parser.add_argument('--verbose', '-v', help='Turn on verbosity', action='store_true', default=False)
	args=parser.parse_args()

	if args.function in ['finderchart', 'fchart', 'fc']:
		assert args.ra != args.dec != None
		FinderChart(args.ra, args.dec)

	elif args.function in ['csv2jsky', 'c2j']:
		assert args.fname != None
		CSVtoJsky(args.fname)

	elif args.function in ['scheduler', 'sched']:
		Scheduler(args.fname, args.DECpenalty)

	elif args.function in ['tnslist', 'tns']:
		MakeTNSlist(args.output, args.dec_limit, args.mag_lower, args.mag_upper, args.time, args.verbose)

	else:
		raise ValueError('Unknown function argument %s' % args.function)