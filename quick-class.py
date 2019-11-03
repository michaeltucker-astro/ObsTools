"""
Michael Tucker 2018
script to help get general classification of 
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.convolution import convolve, Box1DKernel
from matplotlib.widgets import TextBox,CheckButtons,Slider
import ipdb


#Various lines to look for in transient spectra

LINES = {
	'H':[4102.,4340.,4861.,6563.],
	'HeI':[3888.6, 4471.5, 5875.6, 6678.2],
	'HeII':[4686.0],
	'SiII':[6355.0, 5970.,],
	'NaI':[5890.0],
	'FeII':[4325.,4600.,5020.,5200.],
	'SII':[5630.,5450.],
	'CaII':[3950., 8500.],
	'OIII':[4959.0,5007.0,],
	'NII':[6583.0],
	'NeIII':[3869.0],
	'NeV':[3425.0]

}

COLORS = {
	'H':'b',
	'HeI':'g',
	'SiII':'m',
	'NaI':'c',
	'FeII':'m',
	'SII':'c',
	'CaII':'g',
	'HeII':'g',
	'OIII':'r',
	'NII':'g',
	'NeIII':'skyblue',
	'NeV':'skyblue'
}

#Line IDs from the following:
#Ia (Mazzali+ 1993)
#II 
#
TRANSIENTS = {
	'Ia':['FeII','SiII', 'NaI', 'SII','CaII'],
	'Ib/c':[],
	'II':[],
	'CV':['H', 'He'],
	'Mstar':[],
	'AGN':['H', 'He', 'NII', 'NeIII', 'NeV'],

}

#main function call
def main(fname, z, lines, kernel):
	#read SNIFS spectrum 
	wl, fl, err = np.genfromtxt(fname, unpack=True, dtype=float, comments='#')

	if z > 0.5:
		print('Input redshift > 0.5, setting z = 0.5')
		z = 0.5
	if kernel > 50:
		print('Input kernel > 50, setting kernel = 50')
	plot = Plotter(wl, fl, err, kernel, z, lines)
	plot.startPlot()


class Plotter():
	"""
	Interactive plotting class 

	Inputs:
		wl: 1-D wavelength array
		fl: 1-D flux array
		err: 1-D error array
		kernel1: size of smoothing kernel (pixels) for smaller
	"""

	def __init__(self, wl, fl, err, kernel=5, redshift=0.0, lines=[]):
		#assign attributes
		self.wl = wl
		self.fl = fl
		self.err = err
	
		#apply two smoothing kernels
		self.kernel = kernel
		self.smooth = convolve(self.fl, Box1DKernel(self.kernel))

		self.redshift = redshift
		self.lines = [line.strip() for line in lines.split(',')]

	def startPlot(self):
		#setup plot
		self.fig, self.ax = plt.subplots(figsize=(12,8))
		self.ax.set_xlabel('Wavelength (AA)')
		self.ax.set_ylabel('Flux (erg/s/cm2/AA)')

		#actually plot spectrum
		self.spectrum, = self.ax.plot(self.wl, self.fl, 'k-', drawstyle='steps-mid')
		self.plotErrorbar(start=True)
		self.plotSmooth(start=True)
		self.xlims = self.ax.get_xlim()

		self.addLines(self.lines,start=True)
		
		axbox = plt.axes([0.1, 0.1, 0.5, 0.04])
		zslider = Slider(axbox, 'Redshift', 0., 0.5, valinit=self.redshift, valstep=0.001, valfmt='%1.3f')
		zslider.on_changed(self.addRedshift)

		axbox = plt.axes([0.1, 0.06, 0.5, 0.04])
		text_box2 = TextBox(axbox, 'kernel', initial=str(self.kernel))
		text_box2.on_submit(self.resizeKernel)

		axbox = plt.axes([0.1, 0.02, 0.5, 0.04])
		text_box3 = TextBox(axbox, 'lines (csv)', initial=','.join(self.lines))
		text_box3.on_submit(self.addLines)

		rax = plt.axes([0.7, 0.05, 0.14, 0.1])
		labels = ['Obs spec', 'Err spec', 'Smooth spec']
		visibility = [True, False, False]
		check = CheckButtons(rax, labels, visibility)
		check.on_clicked(self.button_press)

		self.ax.set_xlim(self.xlims)
		plt.subplots_adjust(bottom=0.22, top=0.95, left=0.1, right=0.95)

		plt.show()

	def button_press(self, button):
		#what to do when buttons are pushed
		if button == 'Err spec': self.plotErrorbar()
		elif button == 'Obs spec': self.plotSpec()
		elif button == 'Smooth spec': self.plotSmooth()

		self.fig.canvas.draw()

	def plotErrorbar(self, start=False):
		if start:
			upper = self.fl + self.err
			lower = self.fl - self.err
			self.errorbar = self.ax.fill_between(self.wl, lower, upper, color='k', alpha=0.2)
			self.errorbar.set_visible(False)
		else:
			if self.errorbar.get_visible() == False:
				self.errorbar.set_visible(True)
			else:
				self.errorbar.set_visible(False)

	def addLines(self, lines, start=True):
		if start:
			self.SpecLines = {}
			for line in LINES.keys():
				waves = LINES[line]
				allLines = [self.ax.axvline(wave*(1.0+self.redshift), 0.0,1.0,ls='-', lw=2.0, color=COLORS[line]) for wave in waves]
				for ll,wv in zip(allLines, waves):
					ll.set_visible(False)
					ll.restwave = wv

				self.SpecLines[line] = allLines
			if len(lines) == 0: return
		
		if type(self.lines) == str:
			self.lines = [line.strip() for line in lines.split(',')]
		for line in self.lines:
			if line in TRANSIENTS.keys():
				self.lines += TRANSIENTS[line]
				self.lines.remove(line)
			elif line not in LINES.keys(): print('WARNING (addLines): Line %s not found, skipping...')

		for line in LINES.keys():
			if line in self.lines: visible = True
			else: visible = False

			for wave in self.SpecLines[line]: wave.set_visible(visible)

		self.fig.canvas.draw()

	def updateLines(self):
		for line in self.SpecLines.keys():
			for wave in self.SpecLines[line]:
				wave.set_xdata(wave.restwave*(1.0+self.redshift))
		self.fig.canvas.draw()

	def plotSmooth(self, start=False):
		if start:
			self.smooth_plot, = self.ax.plot(self.wl, self.smooth, 'r-')
			self.smooth_plot.set_visible(False)
			return

		visible = self.smooth_plot.get_visible()
		if visible: self.smooth_plot.set_visible(False)
		else: self.smooth_plot.set_visible(True)


	def addRedshift(self, text):
		try:
			self.redshift = float(text)
			self.updateLines()
		except ValueError:
			print('WARNING (addRedshift): Could not convert string to float, check input.')

	def plotSpec(self):
		if self.spectrum.get_visible(): self.spectrum.set_visible(False)
		else: self.spectrum.set_visible(True)

	def resizeKernel(self, text):
		try:
			self.kernel = float(text)
		except ValueError:
			print('WARNING (resizeKernel): Could not convert string to float, check input.')
			return
		self.smooth = convolve(self.fl, Box1DKernel(self.kernel))
		self.smooth_plot.set_ydata(self.smooth)
		self.fig.canvas.draw()


if __name__=='__main__':
	parser = argparse.ArgumentParser(description="quick-class: quick script for guessing classification/redshift of SNIFS transient spex")
	parser.add_argument('fname', help='SNIFS spectrum file', type=str)
	parser.add_argument('--redshift', '-z', help='inital redshift estimate. Default: 0.0', default=0.0, type=float)
	parser.add_argument('--lines', '-l', help='initial lines to plot, csv w/o spaces (e.g. H,He,Si) or a transient type (e.g. Ia/II). Default: None', default='', type=str)
	parser.add_argument('--kernel', '-k', help='pixel width of smoothing kernel. Default: 5', default=5, type=int)

	args = parser.parse_args()
	main(args.fname, args.redshift, args.lines, args.kernel)
