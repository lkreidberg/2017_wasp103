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
        return np.percentile(x, 100. * q)

path = "/home/kreidberg/Projects/Data_reduction/Spitzer/WASP103_11099/Ch2/analysis/run/fgc/ap2750715/fit_phase_curve/2017-01-06_15:28-physical_model/"

files  = [path + "bestfit.pic"]

phasebins = np.linspace(0., 1., 50)
#phasebins = np.linspace(0., 1., 10)
nbins = len(phasebins) - 1
dilution = 1.1587

for f in files:
	gs = gridspec.GridSpec(2, 1, height_ratios = [2.5,1], hspace=0.05)
	ax = plt.subplot(gs[0, 0])


	mcmc = np.load(path + "d-WA103bo21-allparams-trq1lnphysicalbli.npy")
	mcmc = mcmc.T
	mcmc = mcmc[::10000]

	rp, T_s, xi, T_n, delta_T = mcmc[:,1], mcmc[:, 21], mcmc[:, 24], mcmc[:, 25], mcmc[:, 26]

	print "rp, xi, T_n, delta_T", np.median(rp), np.median(xi), np.median(T_n), np.median(delta_T)

	per =  0.925545613
	t0 = 2457080.64041702 
	t = np.linspace(t0, t0+ per, 1000)
	phase = (t-t0)/per - np.floor((t-t0)/per)		

	models = np.zeros((len(t), mcmc.shape[0]))
	for i in range(len(mcmc)): models[:, i] = spiderman_lc.lc(t, 0.115, T_s[i], 4.e-6, 5.e-6, xi[i], T_n[i], delta_T[i], True)
	print "hard coding in rp for spiderman model"

	#bestrp, bestTs, bestxi, bestTn, bestdeltaT = 0.12189, 6.11e3, 0.1446, 667., 2632.
	#plt.plot(t/per, (np.array(spiderman_lc.lc(t, bestrp, bestTs, 4., 5., bestxi, bestTn, bestdeltaT))-1)*1e3, color='purple', linestyle='dotted')
	
	plt.fill_between(phase, (np.apply_along_axis(quantile, 0, models.T, 0.16)-1.)*1e3, (np.apply_along_axis(quantile, 0, models.T, 0.84) - 1.)*1e3, linewidth=0., alpha=0.5, color='r')
	#plt.fill_between(t/per, (np.apply_along_axis(quantile, 0, models.T, 0.025)-1.)*1e3, (np.apply_along_axis(quantile, 0, models.T, 0.975) - 1.)*1e3, linewidth=0., alpha=0.5, color='b')

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


	plt.errorbar(bin_average, (bindata-1.)*1e3*dilution, yerr = binsigma*1e3, fmt = '.k')
	print "phase, err", bin_average[0:5], binsigma[0:5]

	plt.plot(phase, (bestfit-1.)*1e3*dilution, color = 'k', label = "best fit")

	d = ascii.read("../PhaseCurves/PCBands-TiO-NoClouds.dat")
	plt.plot(d['Phase']/360.+0.5, d['Spitzer2']*1.e3, linestyle='dashed', color = '0.5', label="GCM TiO-NoClouds")
	plt.legend()

	#plt.axhline((4.063 + 4.247)/2.*dilution, color='orange')	#eclipses only
#	plt.axhline((4.83 + 4.77)/2.*dilution, color='orange')	#sine curve model

	ax.xaxis.set_major_locator(FixedLocator(np.array([0.1, 0.3, 0.5, 0.7, 0.9])))
	ax.xaxis.set_minor_locator(FixedLocator(np.array([0., 0.2, 0.4, 0.6, 0.8, 1.])))
	ax.set_xticklabels([])
	plt.ylim(-1, 7.)
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

	#plt.savefig("pc_ch2.ps")
	plt.savefig("pc_ch2.png")
