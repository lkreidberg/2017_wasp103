import matplotlib.pyplot as plt
import numpy as np
import pickle
from astropy.io import ascii
import matplotlib.gridspec as gridspec
from pylab import *
from matplotlib import rc
import spiderman_spherical
import transit_lc
from lc_fit import Model, LightCurveData
import spiderman

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def quantile(x, q):
        return np.percentile(x, 100. * q)

fig = plt.figure(figsize = (4.5,4.5))

gs = gridspec.GridSpec(3, 1, hspace=0.1, wspace = 0.05)

# plots HST white light phase curve
bestfits = ["./WFC3_best_fits/bestfit_transit.pic"]

phasebins = np.linspace(-0.15, 0.15, 30)
nbins = len(phasebins) - 1

every_N = 10000

for ii, f in enumerate(bestfits):
	ax = plt.subplot(gs[ii, 0])

	p = pickle.load(open(f, 'rb'))
	d, m, par = p[0], p[1], p[2]		#stores data,  model, and best fit parameters into d, m, & par
	
	ind = d.err < 9.0e7			#indices of outliers
	
	err = d.err[ind]/d.flux[ind]		#normalized per point uncertainty 
	phase = m.phase[ind]
	data_corr = m.data_corr[ind]
#	bestfit = m.bestfit[ind]
	bestfit = m.transit_model[ind]

        ind = phase > 0.5
        phase[ind] -= 1.

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

	plt.errorbar(bin_average, bindata, yerr = binsigma, fmt = '.k', linewidth = 1.)
	plt.plot(phase, bestfit, color = 'b', label = 'best fit', zorder=-1, alpha = 0.8, linewidth = 2.0)

	#formatting

	ax.xaxis.set_major_locator(FixedLocator(np.array([-0.1, 0., 0.1])))
	ax.xaxis.set_minor_locator(FixedLocator(np.array([-0.15, -0.05, 0.05, 0.15])))
	ax.set_xticklabels([])

	#ax.yaxis.set_major_locator(FixedLocator(np.array([0., 0.5, 1, 1.5])))

	plt.ylim(0.985, 1.005)
	plt.xlim(-0.15, 0.15)

	ax.text(0.07, 0.988, 'HST/WFC3\n1.1 - 1.7 $\mu$m', fontsize=10)


#plot Spitzer phase curves

bestfits = ["./Ch1_best_fits/2018-02-07_transit/bestfit.pic", "./Ch2_best_fits/2018_01-31-18_transit/bestfit.pic"]

colors = ['orange', 'red']

per =  0.92555
phasebins = np.linspace(-0.15, 0.15, 30)
nbins = len(phasebins) - 1

for ii, f in enumerate(bestfits):
	ax = plt.subplot(gs[ii+1, 0])

	per =  0.92555
	t0 = 2457080.64041702 
	#t = np.linspace(t0 - 0.1*per, t0+ 0.1*per, 100)
	#phase = (t-t0)/per - np.floor((t-t0)/per)		

	p = pickle.load(open(f, 'rb'))
	data_corr = p[1] 
	err = p[2] 
	bestfit = p[3]
	phase = p[0] 
        phase -= 1.

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

	#print "phase, err", bin_average[0:5], binsigma[0:5]

	plt.errorbar(bin_average, bindata, yerr = binsigma, fmt = '.k', linewidth = 1.)
	plt.plot(phase, bestfit, color = colors[ii], label = 'best fit', alpha = 0.8, linewidth = 2.0)
	#plt.errorbar(phase, data_corr, yerr = err, fmt = '.k', linewidth = 1.)
	#plt.plot(phase, bestfit, color = colors[ii], label = 'best fit', alpha = 0.8, linewidth = 2.0)


	ax.set_xticklabels([])
	plt.ylim(0.985, 1.005)
	plt.xlim(-0.15, 0.15)
        
	
	if ii ==0: 
		plt.ylabel("Relative flux", fontsize = 12, labelpad = 12)
                ax.xaxis.set_major_locator(FixedLocator(np.array([-0.1, 0., 0.1])))
                ax.xaxis.set_minor_locator(FixedLocator(np.array([-0.15, -0.05, 0.05, 0.15])))
		ax.set_xticklabels([])
		ax.text(0.07, 0.988, 'Spitzer Ch 1\n3.6 $\mu$m', fontsize=10)
	if ii == 1: 
		plt.xlabel("Orbital phase", fontsize=12)
                ax.xaxis.set_major_locator(FixedLocator(np.array([-0.1, 0., 0.1])))
                ax.xaxis.set_minor_locator(FixedLocator(np.array([-0.15, -0.05, 0.05, 0.15])))
		ax.set_xticklabels([-0.1, 0.0, 0.1])
		ax.text(0.07, 0.988, 'Spitzer Ch 2\n4.5 $\mu$m', fontsize=10)


gs.tight_layout(fig)
plt.savefig("fig9.pdf")
plt.show()
