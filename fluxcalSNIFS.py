import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
import spectres
from scipy.interpolate import interp1d
import ipdb

def fluxcalSNIFS(spex, plot=True, coord_match=False, verbose=True):

	if verbose: print('\n\nStarting flux calibration routine...')

	names = [spec.object for spec in spex]
	RA = [spec.ra for spec in spex]
	DEC = [spec.dec for spec in spex]
	airmasses = [spec.am for spec in spex]

	STDSidx, STDname, STDwl, STDfl = FindStandards(names, RA, DEC, airmasses, coord_match, verbose)
	if len(STDSidx) == 0: 
		print('\t!!! WARNING (fluxcalSNIFS) !!! No standard stars found, flux calibrations will not be applied!')

	elif len(STDSidx) < 3:
		print('\t WARNING (fluxcalSNIFS): Only %d standard stars found, flux calibration is likely poor.' % len(STDSidx))

	std_am = []
	std_resp = []
	print(STDSidx, STDwl, STDfl)
	for ii, stdwl, stdfl in zip(STDSidx, STDwl, STDfl):
		print(ii)
		spec = spex[ii]
		spec.is_std = True
		resp = GetStdResponse(spec, stdwl, stdfl)
		
		spec.resp = resp.copy()
		spec.fl *= resp
		spec.fluxcal = True

		std_am.append(spec.am)
		std_resp.append(resp)
	
	print(std_am)
	print(std_resp)
	respFunction = MakeRespFunction(std_am, spex[0].wl, std_resp)

	for ii, name, spec, am in zip(range(len(names)), names, spex, airmasses):
		if ii in STDSidx or spec.fluxcal: continue
		else: spec.is_std = False

		if spec.fl.mean() > 1e-5: spec.fl *= 1e-16
		spec.resp_f = respFunction(spec.wl, spec.am)
		ipdb.set_trace()
		spec.fl *= spec.resp_f
		spec.fluxcal = True
		plt.plot(spec.wl, spec.fl/spec.resp_f, 'k-')
		plt.plot(spec.wl, spec.fl, 'r--')
		plt.show()



def FindStandards(names, obsRA, obsDEC, airm, coord_match=False, verbose=True, sep_cut=1.0*3600.0):

	stdFile = '/home/michael/Documents/Science/scripts/Observing/standard_data.dat'
	stdnames, stdRA, stdDEC, stdSpecFiles = np.genfromtxt(stdFile, unpack=True, dtype=str, comments='#')

	if coord_match:
		stdcoords = SkyCoord(stdRA.astype(float)*u.deg, stdDEC.astype(float)*u.deg)
		obscoords = SkyCoord(obsRA*u.deg, obsDEC*u.deg)
	
		if verbose: print('\nMatching coordinates for standard stars...')
		idx, d2d, d3d = match_coordinates_sky(obscoords, stdcoords)
		idx[d2d.arcsec > sep_cut] = np.nan
		if verbose: print('\t%d matches within %3.1f arcsec cut' % (len(d2d[d2d.arcsec < sep_cut]), sep_cut))

	else:
		stdnames = [std.upper() for std in stdnames]
		idx = []
		for name in names:
			if name.upper() in stdnames:
				ii = stdnames.index(name.upper())
				if verbose: print('\t Name-matched %s' % name)
				idx.append(ii)
			else:
				idx.append(np.nan)

	STDwl, STDfl, STDname = [],[],[]
	IDX = np.where(np.isfinite(idx))[0].tolist()
	for ii in idx:
		if ~np.isfinite(ii): continue
		wl, fl = np.genfromtxt(stdSpecFiles[ii], dtype=float, comments='#', unpack=True, usecols=(0,1))
		#spacing = [wl[ii+1] - wl[ii] for ii in range(len(wl)-1)]
		
		#spacing.insert(0, spacing[0])
		#spacing = np.array(spacing)
		STDwl.append(wl)
		#fl /= spacing
		if fl.mean() > 1e-5:
			fl *= 1e-16
		STDfl.append(fl) #convert to erg/cm2/s/A
		STDname.append(stdnames[ii])

	return IDX, STDname, STDwl, STDfl


def GetStdResponse(spec, stdwl, stdfl, spline_order=2, plot=True):
	print(spec.object)
	stdfl_interp = spectres.spectres(spec.wl, stdwl, stdfl)
	resp = stdfl_interp/spec.fl
	default_wl = [3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0, 6500.0, 7100.0, 7500.0, 8000.0, 8500.0, 9000.0, 9500.0]
	spline_wl = [w for w in default_wl if w > spec.wl.min() and w < spec.wl.max()]
	spline_fl = np.interp(spline_wl, spec.wl, resp)

	respFunc = interp1d(spline_wl, spline_fl, kind=spline_order, bounds_error=False, fill_value='extrapolate')
	if plot:
		plt.figure(10)
		plt.plot(spec.wl, spec.fl, 'k-', drawstyle='steps-mid', label='Obs Std Spec')
		plt.plot(stdwl, stdfl, 'r-', alpha=0.7, label='Std Spec Ref')
		plt.plot(spec.wl, spec.fl*respFunc(spec.wl), 'g--', label='Flux-cal obs spec')
		plt.title('Observed and ref standard spectra')
		plt.legend()

		plt.figure(15)
		plt.plot(spec.wl, resp, 'k-', label='Pure resp func')
		plt.plot(spline_wl, spline_fl, 'bo', label='spline points')
		plt.plot(spec.wl, respFunc(spec.wl), 'b-', alpha=0.5, label='response spline fit')
		plt.title('response functions')

		plt.legend()
		plt.show()

	print(spec.wl.shape)

	return respFunc(spec.wl)

class MakeRespFunction():
	def __init__(self, std_am, std_wl, std_resp):
		self.std_am = std_am
		self.std_resp = std_resp
		self.wl = std_wl.copy()

		self.Nstd = len(std_am)

		if self.Nstd < 1: raise ValueError(std_am)
		elif self.Nstd == 1:
			return
		elif self.Nstd == 2:
			self.respFunc = interp1d(np.array(std_am), np.stack(std_resp).T, kind='linear', bounds_error=False, fill_value='extrapolate')
		else:
			self.respFunc = interp1d(np.array(std_am), np.stack(std_resp).T, kind='quadratic', bounds_error=False, fill_value='extrapolate')

	def __call__(self, wl, am):
		print(self.std_resp[0])
		if self.Nstd == 1: 
			ipdb.set_trace()
			interp_resp = np.interp(wl, self.wl, self.std_resp[0])
			return interp_resp
		else:
			resp_interp_am = self.respFunc(am)
			resp_interp_wl = np.interp(wl, self.wl, resp_interp_am)
			return resp_interp_wl



