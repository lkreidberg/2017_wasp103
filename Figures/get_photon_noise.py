import pickle
import numpy as np
import matplotlib.pyplot as plt
import os, glob
from lc_fit import Model, LightCurveData



path  = "WFC3_best_fits/spec_fits/"		#zhang model
files = glob.glob(os.path.join(path, "bestfit*.pic"))		

#phasebins = np.linspace(0.06, 0.44, 5)
#phasebins2 = np.linspace(0.56, 0.94, 5)


for f in files:
	p = pickle.load(open(f, 'rb'))
	data, m = p[0], p[1]			#stores data and model into d & m
	
	ind = data.vis_num < 2

	rms = 1.0e6*np.sqrt(np.mean((m.resid[ind]/data.flux[ind])**2))
        rms_predicted = 1.0e6*np.sqrt(np.mean((data.err[ind]/data.flux[ind])**2))

	print data.wavelength, rms_predicted, rms, rms/rms_predicted
