#import matplotlib
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

def quantile(x, q):
        #return np.percentile(x, [100. * qi for qi in q])
        return np.percentile(x, 100. * q)

#files = ["../PHYSICAL_MODEL_MCMC/ch2_bestfit.pic"]
path = "/home/kreidberg/Projects/Data_reduction/Spitzer/WASP103_11099/Ch1/analysis/run/fgc/ap2500715/fit_phase_curve/2017-01-08_17:01-physical_model/"
files = ["/home/kreidberg/Projects/Data_reduction/Spitzer/WASP103_11099/Ch1/analysis/run/fgc/ap2500715/fit_phase_curve/2017-01-08_17:01-physical_model/bestfit.pic"]

phasebins = np.linspace(0., 1., 50)
nbins = len(phasebins) - 1
dilution = 1.1712

for f in files:
	gs = gridspec.GridSpec(2, 1, height_ratios = [2.5,1], hspace=0.05)
	ax = plt.subplot(gs[0, 0])

	mcmc = np.load(path + "d-WA103bo11-allparams-trq1lnphysicalbliRN.npy").T[::1000]
	#mcmc = np.load(path + "d-WA103bo11-allparams-trq1lnphysicalbli.npy").T[::10000]
	rp, T_s, xi, T_n, delta_T = mcmc[:,1], mcmc[:, 21], mcmc[:, 24], mcmc[:, 25], mcmc[:, 26]

	print "rp, xi, T_n, delta_T", np.median(rp), np.median(xi), np.median(T_n), np.median(delta_T)

	per =  0.92555
	t0 = 2457080.64041702 
	t = np.linspace(t0, t0+ per, 1000)
	phase = (t-t0)/per - np.floor((t-t0)/per)		

	models = np.zeros((len(t), mcmc.shape[0]))
	for i in range(len(mcmc)): models[:, i] = spiderman_lc.lc(t, 0.115, T_s[i], 3.15e-6, 3.95e-6, xi[i], T_n[i], delta_T[i], True)
	print "hard coding in rp for spiderman model"

	
	plt.fill_between(phase, (np.apply_along_axis(quantile, 0, models.T, 0.16)-1.)*1e3, (np.apply_along_axis(quantile, 0, models.T, 0.84) - 1.)*1e3, linewidth=0., alpha=0.5, color='orange')
	#plt.fill_between(t/per, (np.apply_along_axis(quantile, 0, models.T, 0.025)-1.)*1e3, (np.apply_along_axis(quantile, 0, models.T, 0.975) - 1.)*1e3, linewidth=0., alpha=0.5, color='orange')

	p = pickle.load(open(f, 'rb'))
	data_corr = p[1] 
	err = p[2] 
	bestfit = p[3]
	phase = p[0] 
	ind = phase > 1.0
	phase[ind] -= 1.0
	bindata = np.zeros(nbins)
	binsigma = np.zeros(nbins)
	binbestfit = np.zeros(nbins)
	bin_average = np.zeros(nbins)
	
	for j in range(1, len(phasebins)):
		ind = (phase >= phasebins[j-1]) & (phase < phasebins[j])
		if sum(ind)>0:
			bindata[j-1] = sum(data_corr[ind]/err[ind]**2)/sum(1/err[ind]**2) 
			binsigma[j-1] = np.sqrt(1/sum(1/err[ind]**2))
			binbestfit[j-1]=  np.mean(bestfit[ind])
			bin_average[j-1] = (phasebins[j-1]+phasebins[j])/2.

	print "phase, err", bin_average[0:5], binsigma[0:5]

	plt.errorbar(bin_average, (bindata-1.)*1e3*dilution, yerr = binsigma*1e3, fmt = '.k')

	plt.plot(phase, (bestfit-1.)*1e3*dilution, color = 'k', label = 'best fit')

	d = ascii.read("../PhaseCurves/PCBands-TiO-NoClouds.dat")
	plt.plot(d['Phase']/360.+0.5, d['Spitzer1']*1.e3, linestyle='dashed', color = '0.5', label="GCM TiO-NoClouds")
	#d = ascii.read("../PhaseCurves/PCBands-NoTiO-NoClouds.dat")
	#plt.plot(d['Phase']/360.+0.5, d['Spitzer1']*1.e3, linestyle='dotted', color = 'red', label="noTiO-NoClouds")
	#d = ascii.read("../PhaseCurves/PCBands-NoTiO-MgSiO3-10microns.dat")
	#plt.plot(d['Phase']/360.+0.5, d['Spitzer1']*1.e3, linestyle='dotdashed', color = 'orange', label="noTiO-MgSiO3-10um")
	#d = ascii.read("../PhaseCurves/PCBands-TiO-NoClouds-Drag2.dat")
	#plt.plot(d['Phase']/360.+0.5, d['Spitzer1']*1.e3, linestyle='dotted', color = 'green', label="TiO-NoClouds-Drag2")

	plt.legend()


	ax.xaxis.set_major_locator(FixedLocator(np.array([0.1, 0.3, 0.5, 0.7, 0.9])))
	ax.xaxis.set_minor_locator(FixedLocator(np.array([0., 0.2, 0.4, 0.6, 0.8, 1.])))
	ax.set_xticklabels([])
	plt.ylim(-1, 6.)
	plt.ylabel("Planet-to-star flux (ppt)")


	#Plot RESIDUALS 	#
	#########################

	ax = plt.subplot(gs[1, 0])
	binresid = bindata - binbestfit
	plt.axhline(0, color='0.5')
	plt.errorbar(bin_average, binresid*1e3, yerr = binsigma*1e3, fmt = '.k')
	plt.ylim(-2,2)
	
	ax.yaxis.set_major_locator(FixedLocator(np.array([-2,-1,0, 1, 2])))
	ax.set_yticklabels([-2, -1, 0, 1, 2])
	
	ax.xaxis.set_major_locator(FixedLocator(np.array([0.1, 0.3, 0.5, 0.7, 0.9])))
	ax.set_xticklabels([0.1, 0.3, 0.5, 0.7, 0.9])
	ax.xaxis.set_minor_locator(FixedLocator(np.array([0., 0.2, 0.4, 0.6, 0.8, 1.])))
	
	plt.ylabel("Residuals (ppt)")
	plt.xlabel("Orbital phase")

	plt.savefig("pc_ch1.png")
