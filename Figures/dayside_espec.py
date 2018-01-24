import numpy as np
import pickle
import matplotlib.pyplot as plt
import os, glob
from lc_fit import Model, LightCurveData
#import read_spectra_from_Vivien as rs
import matplotlib.gridspec as gridspec
from pylab import *
from matplotlib import rc
import spiderman_lc
from matplotlib import ticker
import seaborn as sns
from astropy.convolution import Gaussian1DKernel, convolve
g = Gaussian1DKernel(stddev=0.8)

sns.set_context("paper", font_scale=1.2)
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})


rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def blackbody(l,T): 
	h = 6.62607e-34   #J/s
        c = 2.997925e8       #m/s
        k = 1.38065e-23     #J/K
        return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T))-1))                 #[W/sr/m^3]

def best_fit_bb(w, y, e, rprs):
	Ts = np.linspace(2900., 3100, 100)
	chibest = 10000.
	Tbest = 0.	
	Tstar = 6110.
	w = np.array(w)

	w, y, e, = w[0:10], y[0:10], e[0:10]

	#get stellar spectrum
	star = np.genfromtxt("wasp103_sed_fluxes.out")
	star_bb = np.interp(w, star[:,0], star[:,1])*1.e24/(w*np.pi*4.)
        #Tstar = 6110.
        #star_bb = blackbody(w*1.0e-6, Tstar)
	outliers = 0.
        chis = []

	for T in Ts:
		model = blackbody(w*1.0e-6, T)/star_bb*rprs**2
		chi2 = np.sum((y - model)**2/e**2)
                chis.append(chi2)
		if chi2 < chibest: 
			chibest, Tbest, outliers = chi2, T, (y-model)/e
	waves_hires = np.linspace(1.0, 5.0, 100)
	star_bb_hires = np.interp(waves_hires, star[:,0], star[:,1])*1.e24/(waves_hires*np.pi*4.)
	#print "Best fit blackbody temp, chisq, and outliers: ", Tbest, chibest/(len(e)-1), outliers

        #find 1-sigma confidence
        idx = (np.abs(chis-(chibest+1.))).argmin()
        """plt.clf()
        plt.plot(Ts, chis)
        plt.axhline(chibest)
        plt.axhline(chibest+1)
        plt.show()"""
        onesigma = Tbest - Ts[idx]
	print "Best fit blackbody temp, chi2: ", Tbest,  "+/-", onesigma, chibest/(len(y) - 1.)
	return waves_hires, blackbody(waves_hires*1.0e-6, Tbest)/star_bb_hires*rprs**2


######################################################################################
# PLOT
plt.figure(figsize = (4,3))

#fits from Mike
data_wl, data, data_err,  best_fit_binned,  wl_hi,  y_hi_best, spec_arr, Tarr, P, samples = pickle.load(open("Mike_models/WASP-103b_grid_DAYSIDE_output.pic"))

print "wavelength, residual from best fit:"
for i in range(len(data_wl)):
    print data_wl[i], (data[i] - best_fit_binned[i])/data_err[i]

plt.errorbar(data_wl, data*1e3, yerr = data_err*1.e3, fmt = 'ow', zorder=1000, ecolor = 'k', markeredgecolor = 'k', markeredgewidth = 1.0)

#plot 1-sigma range
y_minus1sig = np.zeros_like(wl_hi)
y_plus1sig = np.zeros_like(wl_hi)
for i in range(len(wl_hi)): y_minus1sig[i] = np.percentile(spec_arr[:,i], 16)
for i in range(len(wl_hi)): y_plus1sig[i] = np.percentile(spec_arr[:,i], 84)

plt.fill_between(wl_hi, y_minus1sig*1e3, y_plus1sig*1e3, color = 'orange', alpha = 0.5, zorder=-11)
plt.plot(wl_hi, y_hi_best*1e3, zorder = -9, label = 'best fit 1-D model')

#get best fit blackbody
rprs = 0.1127
wave_bb, bb = best_fit_bb(data_wl, data, data_err, rprs)
plt.plot(wave_bb, bb*1e3, linestyle = 'dashed', color = '.5', zorder = -10, label = "blackbody")


#plots gcm
GCM = np.genfromtxt("Vivien_models2/SpectralPC-Phi-TiO-NoClouds-Drag3-NEW-OPA-NEW-PT.dat", delimiter = ",")
plt.plot(GCM[:,0], GCM[:,5]*1e3, color = 'r', label = "$\\tau_\mathrm{drag3}$ GCM")


plt.legend(loc = 'lower right', frameon=True, fontsize = 10)

plt.ylabel("Planet-to-star flux (ppt)")
plt.xlabel("Wavelength (microns)")

plt.xlim(1.1, 1.7)
ymin, ymax = 1.2, 2
plt.ylim(ymin, ymax)

#plots Cartier spectrum
"""s = np.genfromtxt("cartier_dayside_spec.txt")
plt.errorbar(s[:,0], s[:,1]/1.e3, s[:,3]/1.e3, fmt = 'xk')

#plots LK spectrum
s = np.genfromtxt("kreidberg13360_dayside_spec.txt")
plt.errorbar(s[:,0], s[:,1], s[:,2], fmt = 'xr')"""

a = plt.axes([.23, .62, .2, .28]) 

plt.errorbar(data_wl, data*1e3, yerr = data_err*1.e3, xerr = 0.5, fmt = 'ow', zorder=1000, ecolor = 'k', markeredgecolor = 'k', markeredgewidth = 1.0)
plt.fill_between(wl_hi, y_minus1sig*1e3, y_plus1sig*1e3, color = 'orange', alpha = 0.5, zorder=-11)
plt.plot(wl_hi, y_hi_best*1e3, zorder = -9, label = 'best fit')
plt.plot(wave_bb, bb*1e3, linestyle = 'dashed', color = '.5', zorder = -10)

plt.xlim(3, 5)
plt.ylim(3.5,7)


print "NOTE:"
print "you are hard coding rprs for calculating best fit bb"
plt.tight_layout()
plt.savefig("dayside_spectrum.pdf")
#plt.show()
