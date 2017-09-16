import numpy as np
import scipy.interpolate as spi
import matplotlib.pyplot as plt
import matplotlib as mpl
import spiderman as sp
import pickle
from lc_fit import Model, LightCurveData

mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.size'] = 14
mpl.rcParams['lines.linewidth'] = 1.4

models = ["spherical", "zhang", "hotspot_t"]
path = "WFC3_best_fits/"
waverange = [1.1e-6, 1.7e-6]                    #WFC3
colors = ['red', 'orange', 'blue'] 
grays = ['0.7', '0.5', '0.6']
labels = ["Spherical harmonics", "Kinematic", "Two temperature"] 

phasebins = np.linspace(0., 1., 50)
nbins = len(phasebins) - 1

plt.figure(figsize = (6,4))
for i, model in enumerate(models):

	p = pickle.load(open(path+"bestfit_"+model+".pic", "rb")) 

	d, m, par = p[0], p[1], p[2]            #stores data,  model, and best fit parameters into d, m, & par

	dilution = d.dilution + 1.

	ind = d.err < 9.0e7                     #indices of outliers

	err = d.err[ind]/d.flux[ind]            #normalized per point uncertainty 
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

	plt.errorbar(bin_average, (bindata-1.)*1e3*dilution, yerr = binsigma*1e3, linestyle='none', marker = '.', color = colors[i], ecolor = 'k', zorder = 10)

        plt.plot(phase, (bestfit-1.)*1e3*dilution, color = colors[i], label = labels[i])


plt.legend(loc = 'upper left', frameon=False, fontsize=12)
plt.xlabel("Orbital phase")
plt.ylabel("Planet-to-star flux (ppt)")
plt.ylim(-0.2, 2)
plt.xlim(0,1)
plt.tight_layout()
plt.savefig("hst_bestfit.pdf") 
#plt.show()
