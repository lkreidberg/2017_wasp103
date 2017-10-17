import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
#from scipy.interpolate import interp2d

def spectrum(phase, band, model):
	if band == "all":
		#f = "../PhaseCurves/SpectralPC-NoTiO-NoClouds.dat"
		f = "../PhaseCurves/SpectralPC-" + model
		#f = "../PhaseCurves/SpectralPC-TiO-MgSiO3-0.1microns.dat"

		d = ascii.read(f, data_start=1, delimiter=',')
		nrow = len(d)
		waves = np.zeros(nrow)

		for i in range(nrow): waves[i] = d[i][0]

		phases = np.arange(-180, 181, 10)/360. 	+ 0.5	#phase angle in degrees (from file header)
		ncol = len(phases)

		j = 0
		while(phase > phases[j]): j += 1
		print phase, phases[j-1], phases[j]
				
		spec = []
		for i in range(nrow): spec.append( np.interp(phase, np.array([phases[j-1], phases[j]]), np.array([d[i][j-1], d[i][j]]) ) )


		return waves, np.array(spec)
	else:
		f = "../PhaseCurves/PCBands-" + model
		d = ascii.read(f, data_start=1, delimiter=',')
		
		phases = d['Phase']/360.+0.5
		fpfs = d[band]

		return np.interp(phase, phases, fpfs)


