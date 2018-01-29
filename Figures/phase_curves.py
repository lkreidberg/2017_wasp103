import matplotlib.pyplot as plt
import numpy as np
import pickle
from astropy.io import ascii
import matplotlib.gridspec as gridspec
from pylab import *
from matplotlib import rc
import spiderman_lc
import transit_lc
from lc_fit import Model, LightCurveData

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def quantile(x, q):
        return np.percentile(x, 100. * q)

plt.figure(figsize = (7.5,4))

gs = gridspec.GridSpec(3, 2, hspace=0.1, wspace = 0.05)

# plots HST white light phase curve
bestfits = ["/Users/lkreidberg/Desktop/Projects/Observations/HST/WASP103_HST_all/PHYSICAL_MODEL_WHITE_MCMC/bestfit_white.pic"]
mcmcs =  ["/Users/lkreidberg/Desktop/Projects/Observations/HST/WASP103_HST_all/PHYSICAL_MODEL_WHITE_MCMC/mcmc_out_white.npy"]

phasebins = np.linspace(0., 1., 30)
nbins = len(phasebins) - 1

every_N = 10000

for ii, f in enumerate(bestfits):
	ax = plt.subplot(gs[ii, 0])

	mcmc = np.load(mcmcs[ii])
	mcmc = mcmc[::every_N]
	rp, T_s, xi, T_n, delta_T = mcmc[:,0], mcmc[:, 18], mcmc[:, 15], mcmc[:, 16], mcmc[:, 17]

	print "rp, xi, T_n, delta_T", np.median(rp), np.median(xi), np.median(T_n), np.median(delta_T)

	per =  0.925545613
	t0 = 2457080.64041702 
	t = np.linspace(t0, t0+ per, 1000)
	models = np.zeros((len(t), mcmc.shape[0]))
	for i in range(len(mcmc)): models[:, i] = spiderman_lc.lc(t, np.median(rp), T_s[i], 1.1e-6, 1.7e-6, xi[i], T_n[i], delta_T[i], 0., True)*transit_lc.lc(t, np.median(rp))
	#print "hard coding in rp and wavelength for spiderman model"

	phase = (t-t0)/per - np.floor((t-t0)/per)		
	

	p = pickle.load(open(f, 'rb'))
	d, m, par = p[0], p[1], p[2]		#stores data,  model, and best fit parameters into d, m, & par
	
	dilution = d.dilution + 1.

	plt.fill_between(phase, (np.apply_along_axis(quantile, 0, models.T, 0.16)-1.)*1e3*dilution, (np.apply_along_axis(quantile, 0, models.T, 0.84) - 1.)*1e3*dilution, linewidth=0., alpha=0.3, color='blue')

	ind = d.err < 9.0e7			#indices of outliers
	
	err = d.err[ind]/d.flux[ind]		#normalized per point uncertainty 
	phase = m.phase[ind]
	data_corr = m.data_corr[ind]
	bestfit = m.bestfit[ind]

	ind = np.argsort(phase)
	err, phase, data_corr, bestfit = err[ind], phase[ind], data_corr[ind], bestfit[ind] #sorts by phase
	
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

	plt.plot(phase, (bestfit-1.)*1e3*dilution, color = 'k', label = 'best fit')

	ax.xaxis.set_major_locator(FixedLocator(np.array([0.1, 0.3, 0.5, 0.7, 0.9])))
	ax.xaxis.set_minor_locator(FixedLocator(np.array([0., 0.2, 0.4, 0.6, 0.8, 1.])))
	ax.set_xticklabels([])

	ax.yaxis.set_major_locator(FixedLocator(np.array([0., 0.5, 1, 1.5])))

	plt.ylim(-0.2, 1.7)
	plt.xlim(0,1)

	ax.text(0.03, 1.0, 'HST/WFC3\n1.1 - 1.7 $\mu$m', fontsize=10)

	#Plot RESIDUALS 	#
	#########################

	ax = plt.subplot(gs[ii, 1])
	resid = data_corr - bestfit 
	plt.axhline(0, color='0.5', zorder=-10)
	plt.errorbar(bin_average, (bindata - binbestfit)*1e3, yerr = binsigma*1e3, fmt = '.k')
	plt.ylim(-0.2,0.2)
	plt.xlim(0,1)
	
	ax.yaxis.set_major_locator(FixedLocator(np.array([-0.1, 0, 0.1])))
	ax.set_yticklabels([-0.1, 0, 0.1])
	
	ax.xaxis.set_major_locator(FixedLocator(np.array([0.1, 0.3, 0.5, 0.7, 0.9])))
	#ax.set_xticklabels([0.1, 0.3, 0.5, 0.7, 0.9])
	ax.set_xticklabels([])
	ax.xaxis.set_minor_locator(FixedLocator(np.array([0., 0.2, 0.4, 0.6, 0.8, 1.])))
	
	ax.yaxis.tick_right()
	plt.xlabel("Orbital phase", fontsize=12)


#plot Spitzer phase curves

bestfits = ["/Users/lkreidberg/Desktop/Projects/Observations/Spitzer/WASP103_11099/Ch1/analysis/run/fgc_old/ap2500715/fit_phase_curve/2017-01-08_17:01-physical_model/bestfit.pic", "/Users/lkreidberg/Desktop/Projects/Observations/Spitzer/WASP103_11099/Ch2/analysis/run/fgc_old/ap2750715/fit_phase_curve/2017-01-06_15:28-physical_model/bestfit.pic"]

mcmc_output = ["/Users/lkreidberg/Desktop/Projects/Observations/Spitzer/WASP103_11099/Ch1/analysis/run/fgc_old/ap2500715/fit_phase_curve/2017-01-08_17:01-physical_model/d-WA103bo11-allparams-trq1lnphysicalbliRN.npy", "/Users/lkreidberg/Desktop/Projects/Observations/Spitzer/WASP103_11099/Ch2/analysis/run/fgc_old/ap2750715/fit_phase_curve/2017-01-06_15:28-physical_model/d-WA103bo21-allparams-trq1lnphysicalbli.npy"]

colors = ['orange', 'red']

phasebins = np.linspace(0., 1., 30)
nbins = len(phasebins) - 1
dilution = np.array([0.1712, 0.1587])
l1 = np.array([3.15e-6, 4.0e-6])
l2 = np.array([3.95e-6, 5.0e-6])


for ii, f in enumerate(bestfits):
	ax = plt.subplot(gs[ii+1, 0])
	mcmc = np.load(mcmc_output[ii]).T[::every_N]
	rp, T_s, xi, T_n, delta_T = mcmc[:,1], mcmc[:, 21], mcmc[:, 24], mcmc[:, 25], mcmc[:, 26]

	print "rp, xi, T_n, delta_T", np.median(rp), np.median(xi), np.median(T_n), np.median(delta_T)

	per =  0.92555
	t0 = 2457080.64041702 
	t = np.linspace(t0, t0+ per, 1000)
	phase = (t-t0)/per - np.floor((t-t0)/per)		

	models = np.zeros((len(t), mcmc.shape[0]))
	for i in range(len(mcmc)): models[:, i] = spiderman_lc.lc(t, np.median(rp), T_s[i], l1[ii], l2[ii], xi[i], T_n[i], delta_T[i], dilution[ii], True)*transit_lc.lc(t, np.median(rp))
	#print "hard coding in rp for spiderman model"

	plt.fill_between(phase, (np.apply_along_axis(quantile, 0, models.T, 0.16)-1.)*1e3*(1+dilution[ii]), (np.apply_along_axis(quantile, 0, models.T, 0.84) - 1.)*1e3*(1+dilution[ii]), linewidth=0., alpha=0.5, color=colors[ii])

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

	plt.errorbar(bin_average, (bindata-1.)*1e3*(1.+dilution[ii]), yerr = binsigma*1e3, fmt = '.k')
	plt.plot(phase, (bestfit-1.)*1e3*(1.+dilution[ii]), color = 'k', label = 'best fit')


	ax.set_xticklabels([])
	plt.ylim(-1, 6.)
	plt.xlim(0,1)
	
	if ii ==0: 
		plt.ylabel("Planet-to-star flux ($\\times10^{-3}$)", fontsize = 12, labelpad = 10)
		ax.set_xticklabels([])
		ax.text(0.03, 3.8, 'Spitzer Ch 1\n3.6 $\mu$m', fontsize=10)
	if ii == 1: 
		plt.xlabel("Orbital phase", fontsize=12)
		ax.xaxis.set_major_locator(FixedLocator(np.array([0.1, 0.3, 0.5, 0.7, 0.9])))
		ax.set_xticklabels([0.1, 0.3, 0.5, 0.7, 0.9])
		ax.xaxis.set_minor_locator(FixedLocator(np.array([0., 0.2, 0.4, 0.6, 0.8, 1.])))
		ax.text(0.03, 3.8, 'Spitzer Ch 2\n4.5 $\mu$m', fontsize=10)


	#Plot RESIDUALS 	#
	#########################

	ax = plt.subplot(gs[ii+1, 1])
	binresid = bindata - binbestfit
	plt.axhline(0, color='0.5', zorder = -10)
	plt.errorbar(bin_average, binresid*1e3, yerr = binsigma*1e3, fmt = '.k')
	plt.ylim(-2,2)
	plt.xlim(0,1)
	
	ax.yaxis.set_major_locator(FixedLocator(np.array([-1,0, 1])))
	ax.set_yticklabels([-1, 0, 1])
	ax.yaxis.tick_right()
	
	ax.xaxis.set_major_locator(FixedLocator(np.array([0.1, 0.3, 0.5, 0.7, 0.9])))
	ax.set_xticklabels([0.1, 0.3, 0.5, 0.7, 0.9])
	ax.xaxis.set_minor_locator(FixedLocator(np.array([0., 0.2, 0.4, 0.6, 0.8, 1.])))
	
	if ii == 0: 
		ax.set_xticklabels([])
		plt.ylabel("Residuals ($\\times10^{-3}$)", fontsize = 12, labelpad = 10)
		ax.yaxis.set_label_position("right")
	if ii == 1: plt.xlabel("Orbital phase", fontsize=12)


plt.savefig("phase_curves.pdf")
