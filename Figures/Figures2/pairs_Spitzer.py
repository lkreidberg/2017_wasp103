#import matplotlib
#matplotlib.use(.ps')
import matplotlib.pyplot as plt
import numpy as np
import triangle
from pylab import *

def quantile(x, q): return np.percentile(x, [100. * qi for qi in q])
files = ["/home/kreidberg/Projects/Data_reduction/Spitzer/WASP103_11099/Ch1/analysis/run/fgc/ap2500715/fit_phase_curve/2017-01-08_17:01-physical_model/d-WA103bo11-allparams-trq1lnphysicalbliRN.npy", "/home/kreidberg/Projects/Data_reduction/Spitzer/WASP103_11099/Ch2/analysis/run/fgc/ap2750715/fit_phase_curve/2017-01-06_15:28-physical_model/d-WA103bo21-allparams-trq1lnphysicalbli.npy"]

outnames = ["ch1pairs.png", "ch2pairs.png"]
	
indices = [[0, 1, 4, 6, 8, 21, 24, 25, 26, 30, 31, 32], [0, 1, 4, 6, 8, 21, 24, 25, 26]]
labels = [["trmid", "rp", "c", "u1", "v", "Ts", "xi", "Tn", "deltaT", "white", "red", "gamma"], ["trmid", "rp", "c", "u1", "v", "Ts", "xi", "Tn", "deltaT"]]

#files = ["/home/kreidberg/Projects/Data_reduction/Spitzer/WASP103_11099/Ch2/analysis/run/fgc/ap2750715/fit_phase_curve/2017-01-02_11:54-physical_model/d-WA103bo21-allparams-trq1lnphysicalbli.npy", "/home/kreidberg/Projects/Data_reduction/Spitzer/WASP103_11099/Ch2/analysis/run/fgc/ap2750715/fit_phase_curve/2017-01-02_16:44-physical_model_w_quadramp/d-WA103bo21-allparams-trq1qdphysicalbli.npy"]

#indices = [[0, 1, 4, 6, 8, 21, 24, 25, 26], [0, 1, 4, 6, 8, 9, 22, 25, 26, 27]]
#labels = [["trmid", "rp", "c", "u1", "v", "Ts", "xi", "Tn", "deltaT"], ["trmid", "rp", "c", "u1", "v", "v1", "Ts", "xi", "Tn", "deltaT"]]

#outnames = ["ch2pairs.png","ch2pairs_quad.png"]

#inds = [[25],[26]]
#inds2 = [[26],[27]]

"""files = ["/home/kreidberg/Projects/Data_reduction/Spitzer/WASP103_11099/Ch2/analysis/run/fgc/ap2750715/fit_phase_curve/2017-01-06_15:28-physical_model/d-WA103bo21-allparams-trq1lnphysicalbli.npy"]
indices = [[0, 1, 4, 6, 8, 21, 24, 25, 26]]
labels = [["trmid", "rp", "c", "u1", "v", "Ts", "xi", "Tn", "deltaT"]]

inds = [[25]]
inds2 = [[26]]
outnames = ["ch2pairs.png"]

for i, f in enumerate(files):
	d = np.load(f)
	Tn = d[inds[i]].T
	deltaT = d[inds2[i]].T

	q = np.array([0.14, 0.5, 0.86])
	qTn, qdeltaT = quantile(Tn, q), quantile(deltaT, q) 
	print qTn[2] - qTn[0], qdeltaT[2] - qdeltaT[0]"""

for i, f in enumerate(files):
	d = np.load(f)
	d = d[indices[i]].T[::100]
	hist(d[:,1])
	plt.show()
	fig = triangle.corner(d, labels = labels[i], range = 0.9, show_titles=True) 
	#plt.tight_layout()
	plt.savefig(outnames[i])
