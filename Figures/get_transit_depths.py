import numpy as np
import pickle
import glob
import os
from lc_fit import Model, LightCurveData


path  = "WFC3_best_fits/spec_fits_transit/"		#zhang model
files = glob.glob(os.path.join(path, "bestfit*.pic"))		

for f in files:
	p = pickle.load(open(f, 'rb'))
	data, m, par = p[0], p[1], p[2]			#stores data and model into d & m
        rps = par[data.par_order['rp']*data.nvisit:(1 + data.par_order['rp'])*data.nvisit]
        rp = rps[0]
        depth = rp**2
        print data.wavelength, depth
