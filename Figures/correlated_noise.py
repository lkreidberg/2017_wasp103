import numpy as np
import matplotlib.pyplot as plt
from lc_fit_nospiderman import Model, LightCurveData
import pickle

path = './WFC3_best_fits/spec_fits/'
bestfits = ['bestfit_1.18.pic']

for ii, f in enumerate(bestfits):
    p = pickle.load(open(path + f, 'rb'))
    d, m, par = p[0], p[1], p[2] 
    
    ind = d.err < 9.0e7                     #indices of outliers

    err = d.err[ind]/d.flux[ind]            #normalized per point uncertainty 
    phase = m.phase[ind]
    data_corr = m.data_corr[ind]
    bestfit = m.bestfit[ind]

    """ind = np.argsort(phase)
    err, phase, data_corr, bestfit = err[ind], phase[ind], data_corr[ind], bestfit[ind] #sorts by phase

    bins = np.arange(10, 100, 2)
    for nbins in bins:
        phasebins = np.linspace(0., 1., nbins)

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

        print nbins, np.mean(binsigma)"""

    nexp = 636 
    nperbin = np.arange(10, 60)

    for n in nperbin:
        nbins = int(nexp/n)
        binsigma = np.zeros(nbins)
        for i in range(nbins):
            binsigma[i] = np.sqrt(1/sum(1/err[i*n:(i+1)*n]**2))

        print n, np.mean(binsigma)



