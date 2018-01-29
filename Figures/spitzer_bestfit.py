import numpy as np
import scipy.interpolate as spi
import matplotlib.pyplot as plt
import matplotlib as mpl
import spiderman as sp
import pickle
#from lc_fit import Model, LightCurveData

mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.size'] = 12
mpl.rcParams['lines.linewidth'] = 1.4


def initialize_params(spider_params):
    spider_params.n_layers= 5
    spider_params.t0= 200               # Central time of PRIMARY transit [days]
    spider_params.per= 0.81347753       # Period [days]
    spider_params.a_abs= 0.01526        # The absolute value of the semi-major axis [AU]
    spider_params.inc= 82.33            # Inclination [degrees]
    spider_params.ecc= 0.0              # Eccentricity
    spider_params.w= 90                 # Argument of periastron
    spider_params.rp= 0.1594            # Planet to star radius ratio
    spider_params.a= 4.855              # Semi-major axis scaled by stellar radius
    spider_params.p_u1= 0               # Planetary limb darkening parameter
    spider_params.p_u2= 0               # Planetary limb darkening parameter
    return spider_params

models = ["spherical", "zhang", "hotspot_t"]
path = "Ch2_best_fits/"
bestfits = ["2017-09-14_18:26-spherical/", "2017-09-15_11:42-zhang/", "2017-09-15_11:45-hotspot/"]
waverange = [4.0e-6, 5.0e-6]                    #Spitzer Ch 2
colors = ['orange', 'red', 'blue'] 
grays = ['0.7', '0.5', '0.6']
labels = ["Spherical harmonics", "Kinematic", "Two temperature"]

phasebins = np.linspace(0., 1., 50)
nbins = len(phasebins) - 1
dilution =  0.1587                  #spitzer ch 2
#dilution = 0.1712                   #spitzer ch 1

#plt.figure(figsize = (6,4))
plt.figure(figsize = (3.5,2.2))
for i, model in enumerate(models):
    spider_params = sp.ModelParams(brightness_model=model)
    spider_params = initialize_params(spider_params)
    spider_params.l1 = waverange[0]
    spider_params.l2 = waverange[1]

    p = pickle.load(open(path+bestfits[i]+"bestfit.pic", "rb")) 

    if model == "zhang":
	spider_params.xi= 0.3       # Ratio of radiative to advective timescale
        spider_params.T_n= 700	# Temperature of nightside
        spider_params.delta_T= 900 # Day-night temperature contrast
        spider_params.T_s = 4500    # Temperature of the star

    elif model == "spherical":
        spider_params.degree=2
        spider_params.la0 = 0
        spider_params.lo0 = 0
        spider_params.sph = [2.75677373e-03, 2.33711161e-04, 1.88579054e-23, 1.21500964e-03]

    elif model == "hotspot_t":
        print "hi"

    else:
        print "model not supported"

    #t= spider_params.t0 + np.linspace(0, + spider_params.per,100)
    #lc = spider_params.lightcurve(t)

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

    if i==1:
        plt.errorbar(bin_average, (bindata-1.)*1e3*(1.+dilution), yerr = binsigma*1e3, linestyle='none', marker = '.', color = colors[i], ecolor = 'k', zorder = 10)
        #plt.plot(phase, (bestfit-1.)*1e3*(1.+dilution), label = labels[i], color = colors[i])
        plt.plot(phase, (bestfit-1.)*1e3*(1.+dilution), label = labels[i], color = 'k') 

#plt.legend(loc = 'upper left', frameon=False, fontsize=12)
plt.xlabel("Orbital phase")
plt.ylabel("Planet-to-star flux ($\\times10^{-3}$)")
#plt.ylim(-0.7, 8)
plt.ylim(-0.7, 6)
plt.xlim(0,1)
#plt.show()
plt.tight_layout()
plt.savefig("spitzer_bestfit.pdf") 
