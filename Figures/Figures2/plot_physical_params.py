import matplotlib
#matplotlib.use('ps')
import matplotlib.pyplot as plt
import numpy as np
import pickle
from astropy.io import ascii
import matplotlib.gridspec as gridspec
from pylab import *
from matplotlib import rc
import spiderman_lc

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
	
mcmc = np.load("/home/kreidberg/Projects/Data_reduction/Spitzer/WASP103_11099/Ch2/analysis/run/fgc/ap2750715/fit_phase_curve/2016-12-24_14:37-physical_model/d-WA103bo21-allparams-trq1lnphysicalbli.npy")
mcmc = mcmc.T
mcmc = mcmc[::100]

rp, T_s, xi, T_n, delta_T = mcmc[:,1], mcmc[:, 21], mcmc[:, 24], mcmc[:, 25], mcmc[:, 26]


bestrp, bestTs, bestxi, bestTn, bestdeltaT = 0.12189, 6.11e3, 0.1446, 667., 2632. 

plt.subplot(221)
hist(rp)
plt.axvline(bestrp)

plt.subplot(222)
hist(T_s)
plt.axvline(bestTs)

plt.subplot(223)
hist(T_n)
plt.axvline(bestTn)

plt.subplot(224)
hist(delta_T)
plt.axvline(bestdeltaT)

plt.show()
