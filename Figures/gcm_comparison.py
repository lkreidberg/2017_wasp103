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
import os, glob


rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def quantile(x, q):
        return np.percentile(x, 100. * q)

plt.figure(figsize = (7.5,4))

gs = gridspec.GridSpec(2, 1, hspace=0.1, wspace = 0.05)

# plots HST white light phase curve
bestfits = ["./WFC3_best_fits/bestfit_spherical.pic"]

gcm_path = "./GCM_From_Vivien/GCMmodels"
gcms = glob.glob(os.path.join(gcm_path, "PC*PT.dat"))

phasebins = np.linspace(0., 1., 30)
nbins = len(phasebins) - 1

every_N = 10000

for ii, f in enumerate(bestfits):
	ax = plt.subplot(gs[0, ii])


	per =  0.925545613
	t0 = 2457080.64041702 
	t = np.linspace(t0, t0+ per, 1000)

	phase = (t-t0)/per - np.floor((t-t0)/per)		
	
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

	plt.errorbar(bin_average, (bindata-1.)*1e3*dilution, yerr = binsigma*1e3, fmt = '.k')

	plt.plot(phase, (bestfit-1.)*1e3*dilution, color = 'k', label = 'best fit')

	ax.xaxis.set_major_locator(FixedLocator(np.array([0.1, 0.3, 0.5, 0.7, 0.9])))
	ax.xaxis.set_minor_locator(FixedLocator(np.array([0., 0.2, 0.4, 0.6, 0.8, 1.])))
	ax.set_xticklabels([])

	ax.yaxis.set_major_locator(FixedLocator(np.array([0., 0.5, 1, 1.5])))

	plt.ylim(-0.2, 1.7)
	plt.xlim(0,1)

	ax.text(0.03, 1.0, 'HST/WFC3\n1.1 - 1.7 $\mu$m', fontsize=10)


        #plot gcms
        for g in gcms:
            model = np.genfromtxt(g, delimiter = ',') 
            plt.plot(model[:,0]/np.max(model[:,0])/2. + 0.5, model[:,5]*1e3)
            

	#Plot RESIDUALS 	#
	#########################

	ax = plt.subplot(gs[1, ii])
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



plt.savefig("gcm_comparison.pdf")
