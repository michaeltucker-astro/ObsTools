import numpy as np
import matplotlib.pyplot as plt
import requests
import os, sys, glob
from bs4 import BeautifulSoup as BS

ASASSN_HTML = 'http://www.astronomy.ohio-state.edu/asassn/transients.html'
NCOL = 12

def main(query='asassn', otype='sn', declim=-40., verb=1):
	if query == 'asassn':
		objects = QueryASASSN(otype=otype, declim=declim, verb=verb)

	

def QueryASASSN(otype='sn', declim=-40., verb=1):
	objects = parseHTMLtable(ASASSN_HTML)
	trim_objs = [obj for obj in objects if obj.otype is None and obj.decdeg >= declim and obj.ctype=='SN' and obj.discdate.split('-')[0] == '2019']
	if verb: print('ASASSN: %d good objs trimmed from %d total objs.' % (len(trim_objs), len(objects)))



class Object():
	def __init__(self, ATname, DiscName, ra, dec, discmag=None, discdate=None, curmag=None, fchart=None, otype=None, ctype=None, hostz=None):
		self.ATname=ATname
		self.name=DiscName
		self.ra = ra
		self.dec = dec
		self.discmag = discmag
		self.discdate=discdate
		self.curmag=curmag
		self.fchart=fchart
		self.otype=otype #object type
		self.ctype=ctype # candidate type
		self.hostz = hostz
		self.convert_coords()

	def convert_coords(self):
		rastr = self.ra.split(':')
		self.radeg = float(rastr[0])*15. + float(rastr[1])/60. + float(rastr[2])/3600.
		if self.dec[0] == '-':
			decstr = self.dec[1:].split(':')
			sign = -1.
		else:
			decstr = self.dec.split(':')
			sign = 1.
		decdeg = float(decstr[0]) + float(decstr[1])/60. + float(decstr[2])/3600.
		self.decdeg = decdeg*sign


def parseHTMLtable(url):
	resp = requests.get(url)
	soup = BS(resp.content.decode(), 'lxml')
	table = soup.find_all('table')[0]
	table_data = table.find_all('td')[4:]


	objects = []
	for ii in range(0, len(table_data), NCOL):
		row = table_data[ii:ii+NCOL]
		obj = makeObjectFromRow(row)
		objects.append(obj)

	return objects

def makeObjectFromRow(row):
	name, othername, TNSurl, ra, dec, date, mag, SDSSurl, DSSurl, VIZurl, otype, comments = [np.str.replace(str(ele)[4:-5], ',', ' ').strip() for ele in row]
	if 'SN candidate' in comments:
		ctype='SN'
		if 'Type unknown' in comments: otype=None
		else: 
			com_data = comments.split()
			try:
				idx = com_data.index('Type')
				otype = com_data[idx+1]
			except:
				otype=None
		if 'z unknown' in comments:
			hostz = None
		elif 'hostz=' in comments:
			com_data = comments.split()
			hostz_str = [item.strip() for item in com_data if 'hostz=' in item][0]
			hostz = (hostz_str.split('=')[-1].strip())
		elif 'z=' in comments:
			com_data = comments.split()
			hostz_str = [item.strip() for item in com_data if 'z=' in item][0]
			hostz = (hostz_str.split('=')[-1].strip())
		else:
			hostz=None

	elif 'CV candidate' in comments:
		ctype='CV'
		otype=None
		hostz=None
	elif 'microlensing' in comments.lower():
		ctype='microlensing'
		otype=None
		hostz=None
	else:
		ctype=None
		otype=None
		hostz=None
	
	try:
		float(mag)
	except:
		mag=0.0

	if hostz is not None:
		try:
			hostz = float(hostz)
		except:
			hostz = float(hostz[:-1])

	obj = Object(name, othername, ra, dec, float(mag), date, None, None, otype, ctype, hostz)
	return obj


if __name__=='__main__':
	main()
