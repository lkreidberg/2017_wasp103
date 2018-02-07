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
files = glob.glob(os.path.join(path, "bestfit*.pic"))		
for f in files:
	p = pickle.load(open(f, 'rb'))
	d, m = p[0], p[1]			

	x, y = m.phase, m.bestfit_no_eclipse
	ind = np.argsort(x)
	x, y = x[ind], y[ind]
	

        #ind = m.bestfit_no_eclipse == np.max(m.bestfit_no_eclipse)


	fnew = interp1d(x, y, kind = 'cubic')
	xnew = np.linspace(0.4, 0.6, 1000)
	ynew = fnew(xnew)

	plt.plot(xnew, ynew, label = f)
	ind = ynew == np.max(ynew)
	offset = xnew[ind]
        print f, (offset - 0.5)*180./np.pi

plt.legend()
plt.show()
