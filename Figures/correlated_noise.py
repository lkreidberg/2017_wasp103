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

    nexp = len(phase)
    """binsz = np.arange(5, 100, 1)
    rms = np.zeros_like(binsz)

    for k, n in enumerate(binsz):
        nbins = int(nexp/n)
        resid = np.zeros(nbins)
        for i in range(nbins):
            ind1, ind2 = i*n, (i+1)*n
            resid[i] = np.mean(data_corr[ind1:ind2]) - np.mean(bestfit[ind1:ind2])

        rms[k] = np.sqrt(np.mean((resid)**2))*1e6

    plt.plot(binsz, rms, color = 'r')
    nbins = nexp/binsz
    #plt.errorbar(binsz, rms, rms/np.sqrt(2*nbins), fmt = '.k')"""

    nbins = np.arange(5, 100, 10)
    rms = np.zeros_like(nbins)

    for k, nbin in enumerate(nbins):
        resid = np.zeros(nbin)
        for i in range(nbin):
            ppb = int(nexp/nbin)
            ind1, ind2 = i*ppb, (i+1)*ppb
            resid[i] = np.mean(data_corr[ind1:ind2]) - np.mean(bestfit[ind1:ind2])
            
        rms[k] = np.sqrt(np.mean((resid)**2))*1e6

    plt.plot(nbins, rms, color = 'r')

    plt.gca().set_xscale('log')
    plt.gca().set_yscale('log')

plt.show()
