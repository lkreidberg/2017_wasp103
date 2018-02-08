import pickle
from astropy.io import ascii
import transit_lc
from lc_fit import Model, LightCurveData
import os, glob
import numpy as np


# plots HST white light phase curve
#bestfits = ["./WFC3_best_fits/old_white_fits/bestfit_zhang_allsys.pic"]
bestfits = ["./WFC3_best_fits/bestfit_spherical.pic"]

for ii, f in enumerate(bestfits):
	p = pickle.load(open(f, 'rb'))
	d, m, par = p[0], p[1], p[2]		#stores data,  model, and best fit parameters into d, m, & par
	

	ind = d.err < 9.0e7			#indices of outliers
	
	err = d.err[ind]/d.flux[ind]		#normalized per point uncertainty 
	phase = m.phase[ind]
	data_corr = m.data_corr[ind]
	allsys = m.lc[ind]	
	t = d.time[ind]
	visit_sys = m.all_sys[ind]

	print "WFC3 white: obs, exp rms (ppm)", np.std(m.resid[ind]/d.flux[ind])*1e6, np.sqrt(1./np.median(d.flux[ind]))*1e6


fluxconv = [306.126, 266.648]
#calculated from Jonathan Fraine's code https://github.com/exowanderer/ExoplanetTSO/blob/master/ExoplanetTSO_Auxiliary.py
"""fluxConv  = testheader['FLUXCONV']
expTime   = testheader['EXPTIME']
gain      = testheader['GAIN']
fluxConversion = expTime*gain / fluxConv"""

bestfits = ["Ch1_best_fits/2018-02-07_14:28-spherical/bestfit.pic",  "Ch2_best_fits/2018-02-07_14:07-spherical/bestfit.pic"]

for ii, f in enumerate(bestfits):
	p = pickle.load(open(f, 'rb'))
	data_corr = p[1] 
	err = p[2] 
	bestfit = p[3]
	phase = p[0] 

	abscissauc = p[10]
	binfluxuc = p[11]
	binstduc = p[12]
	bestfit = p[13]
	abscissa = p[15]
	sys = p[16]

	print "rms obs, exp (ppm)", 1.0e6*np.std((binfluxuc - bestfit)/binfluxuc), 1.e6/np.sqrt(np.median(binfluxuc*fluxconv[ii]))

#calculates rms for WFC3 spectroscopic best fits
path = "./WFC3_best_fits/spec_fits/"
files = glob.glob(os.path.join(path, "*"))	
for f in files: 
	p = pickle.load(open(f, 'rb'))
        d, m, par = p[0], p[1], p[2]            #stores data,  model, and best fit parameters into d, m, & par

        ind = d.err < 9.0e7                     #indices of outliers

        print "WFC3 spec: obs, exp rms (ppm)", d.wavelength, m.rms, 1.0e6*np.sqrt(np.mean((d.err[ind]/d.flux[ind])**2)),  m.rms/(1.0e6*np.sqrt(np.mean((d.err[ind]/d.flux[ind])**2)))

