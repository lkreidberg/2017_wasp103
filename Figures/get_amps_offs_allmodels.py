import pickle
import numpy as np
import matplotlib.pyplot as plt
import os, glob
from lc_fit import Model, LightCurveData
import matplotlib.gridspec as gridspec
from pylab import *
from astropy.convolution import Gaussian1DKernel, convolve
from scipy.interpolate import interp1d

path  = "WFC3_best_fits/"
models = ["hotspot_t", "zhang", "spherical", "sincos"]
for i, model in enumerate(models):
	p = pickle.load(open(path+"bestfit_"+model+".pic", "rb")) 
	d, m = p[0], p[1]			

	x, y = m.phase, m.bestfit_no_eclipse
	ind = np.argsort(x)
	x, y = x[ind], y[ind]
	

	fnew = interp1d(x, y, kind = 'cubic')
	xnew = np.linspace(0.4, 0.6, 1000)
	ynew = fnew(xnew)

	plt.plot(xnew, ynew, label = f)
	ind = ynew == np.max(ynew)
	offset = xnew[ind]
        print f, (offset - 0.5)*180./np.pi

#plt.legend()
#plt.show()
