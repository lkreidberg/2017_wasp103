import pickle
import numpy as np
import matplotlib.pyplot as plt
import os, glob
from lc_fit import Model, LightCurveData
import matplotlib.gridspec as gridspec
from pylab import *
from astropy.convolution import Gaussian1DKernel, convolve

path  = "WFC3_best_fits/"
files = glob.glob(os.path.join(path, "bestfit*.pic"))		
for f in files:
	p = pickle.load(open(f, 'rb'))
	d, m = p[0], p[1]			
        ind = m.bestfit == np.max(m.bestfit)
        print m.phase[ind]
