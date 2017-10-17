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
from lc_fit import Model, LightCurveData

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def quantile(x, q):
        #return np.percentile(x, [100. * qi for qi in q])
        return np.percentile(x, 100. * q)


files = ["../PHYSICAL_MODEL_WHITE_MCMC/bestfit_white.pic"]

phasebins = np.linspace(0., 1., 40)
nbins = len(phasebins) - 1

for f in files:
	gs = gridspec.GridSpec(2, 1, height_ratios = [2.5,1], hspace=0.05)
	ax = plt.subplot(gs[0, 0])

	mcmc = np.load("../PHYSICAL_MODEL_WHITE_MCMC/mcmc_out_white.npy")
	mcmc = mcmc[::10000]
	rp, T_s, xi, T_n, delta_T = mcmc[:,0], mcmc[:, 18], mcmc[:, 15], mcmc[:, 16], mcmc[:, 17]

	print "rp, xi, T_n, delta_T", np.median(rp), np.median(xi), np.median(T_n), np.median(delta_T)

	per =  0.925545613
	t0 = 2457080.64041702 
	t = np.linspace(t0, t0+ per, 1000)
	models = np.zeros((len(t), mcmc.shape[0]))
	for i in range(len(mcmc)): models[:, i] = spiderman_lc.lc(t, 0.115, T_s[i], 1.1e-6, 1.7e-6, xi[i], T_n[i], delta_T[i], True)
	print "hard coding in rp and wavelength for spiderman model"

	phase = (t-t0)/per - np.floor((t-t0)/per)		
	
	plt.fill_between(phase, (np.apply_along_axis(quantile, 0, models.T, 0.16)-1.)*1e3, (np.apply_along_axis(quantile, 0, models.T, 0.84) - 1.)*1e3, linewidth=0., alpha=0.3, color='blue')
#	plt.fill_between(phase, (np.apply_along_axis(quantile, 0, models.T, 0.025)-1.)*1e3, (np.apply_along_axis(quantile, 0, models.T, 0.975) - 1.)*1e3, linewidth=0., alpha=0.2, color='blue')

	p = pickle.load(open(f, 'rb'))
	d, m, par = p[0], p[1], p[2]		#stores data,  model, and best fit parameters into d, m, & par
	
	dilution = d.dilution + 1.


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

	#plt.errorbar(phase, (data_corr-1.)*1e3*dilution, yerr = err*1e3, fmt = '.k')
	plt.errorbar(bin_average, (bindata-1.)*1e3*dilution, yerr = binsigma*1e3, fmt = '.k')

	plt.plot(phase, (bestfit-1.)*1e3*dilution, color = 'k', label = 'best fit')

	"""nvisit = 4
	best_T_s, best_xi, best_T_n, best_delta_T =  par[d.par_order['T_s']*nvisit], par[d.par_order['xi']*nvisit], par[d.par_order['T_n']*nvisit],  par[d.par_order['delta_T']*nvisit] 
	phase = (t - t0)/per - np.floor((t - t0)/per) 
	plt.plot(phase, (spiderman_lc.lc(t, 0.115, best_T_s, 1.1e-6, 1.7e-6, best_xi, best_T_n, best_delta_T, True) - 1.)*1.e3, color='r')""" #checking to see whether bestfit is right

	"""d = ascii.read("../PhaseCurves/PCBands-TiO-NoClouds.dat")
	plt.plot(d['Phase']/360.+0.5, d['WFC3']*1.e3, linestyle='dashed', color = '0.5', label="GCM TiO-NoClouds")

	d = ascii.read("../PhaseCurves/PCBands-TiO-NoClouds-Drag1.dat")
	plt.plot(d['Phase']/360.+0.5, d['WFC3']*1.e3, linestyle='dotted', color = '0.5', label="Drag 1")

	d = ascii.read("../PhaseCurves/PCBands-TiO-NoClouds-Drag2.dat")
	plt.plot(d['Phase']/360.+0.5, d['WFC3']*1.e3, linestyle='dashdot', color = '0.5', label="Drag 2")
	plt.legend(loc=2, fontsize=9)"""

	ax.xaxis.set_major_locator(FixedLocator(np.array([0.1, 0.3, 0.5, 0.7, 0.9])))
	ax.xaxis.set_minor_locator(FixedLocator(np.array([0., 0.2, 0.4, 0.6, 0.8, 1.])))
	ax.set_xticklabels([])
	#plt.ylim(-0.2, 2.)
	plt.ylim(-0.2, 1.7)
	plt.ylabel("Planet-to-star flux (ppt)")


	#Plot RESIDUALS 	#
	#########################

	ax = plt.subplot(gs[1, 0])
	resid = data_corr - bestfit 
	plt.axhline(0, color='0.5')
	#plt.errorbar(phase, resid*1e3, yerr = err*1e3, fmt = '.k')
	plt.errorbar(bin_average, (bindata - binbestfit)*1e3, yerr = binsigma*1e3, fmt = '.k')
	plt.ylim(-0.2,0.2)
	
	ax.yaxis.set_major_locator(FixedLocator(np.array([-0.1, 0, 0.1])))
	ax.set_yticklabels([-0.1, 0, 0.1])
	
	ax.xaxis.set_major_locator(FixedLocator(np.array([0.1, 0.3, 0.5, 0.7, 0.9])))
	ax.set_xticklabels([0.1, 0.3, 0.5, 0.7, 0.9])
	ax.xaxis.set_minor_locator(FixedLocator(np.array([0., 0.2, 0.4, 0.6, 0.8, 1.])))
	
	plt.ylabel("Residuals (ppt)")
	plt.xlabel("Orbital phase")

	#plt.savefig("pc_hst.ps")
	plt.savefig("pc_hst.png")
