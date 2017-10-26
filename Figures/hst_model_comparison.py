import numpy as np
import scipy.interpolate as spi
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
import spiderman as sp
import pickle
from lc_fit import Model, LightCurveData
import matplotlib.gridspec as gridspec
import seaborn as sns

mpl.rcParams['font.sans-serif'] = 'Arial'

sns.set_context("notebook", font_scale = 1.0)
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"out", "ytick.direction":"out"})


#return blackbody spectral radiance as a function of wavelength lambda (in m) and temperature T (Kelvin)
def bb(l, T):
	h = 6.6261e-34   #J/s
	c = 2.9979e8       #m/s
	k = 1.38065e-23     #J/K
	return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T)-1)))

def get_T(relative):
	l = 1.4e-6
	Ts = 6110.
	rprs = 0.11092
	T = np.arange(10, 4000, 10)
	rel = bb(l, T)/bb(l, Ts)*rprs**2
	diff = np.abs(rel - relative)
	return T[diff == np.min(diff)]

def convert_to_T(flux):
	Ts = 6110. 
	l = 1.4e-6
	n, m = flux.shape
	T = np.zeros_like(flux)
	for i in range(m):
		for j in range(n):
			relative = flux[i,j]
			T[i,j] = get_T(relative)
			#if relative > 0: T[i,j] = get_T(relative)
			#else: 
                         #       T[i,j] = -1000.
	return T

def initialize_params(spider_params):
    spider_params.n_layers= 5
    spider_params.t0= 200               # Central time of PRIMARY transit [days]
    spider_params.per= 0.81347753       # Period [days]
    spider_params.a_abs= 0.01526        # The absolute value of the semi-major axis [AU]
    spider_params.inc= 82.33            # Inclination [deg]
    spider_params.ecc= 0.0              # Eccentricity
    spider_params.w= 90                 # Argument of periastron
    spider_params.rp= 0.111            # Planet to star radius ratio
    spider_params.a= 4.855              # Semi-major axis scaled by stellar radius
    spider_params.p_u1= 0               # Planetary limb darkening parameter
    spider_params.p_u2= 0               # Planetary limb darkening parameter
    spider_params.T_s = 6110.		# Temperature of the star
    return spider_params

models = ["spherical", "zhang", "hotspot_t"]
labels = ["Sph. harmonics", "Kinematic", "Two temp."] 
path = "WFC3_best_fits/"
waverange = [1.1e-6, 1.7e-6]                    #Spitzer Ch 2

fig =plt.figure(figsize = (10,3.33))
gs = gridspec.GridSpec(3, 5, width_ratios = [3, 0.1, 1,1,0.1], height_ratios = [1,1, 0.01]) 
ncol = 2

for i, model in enumerate(models):	
    row, col = int(np.floor(i/ncol)), i%ncol
    col += 2
    ax = plt.subplot(gs[row, col])

    if model == "spherical": spider_params = sp.ModelParams(brightness_model=model, thermal=True)
    else: spider_params = sp.ModelParams(brightness_model=model)

    spider_params = initialize_params(spider_params)
    spider_params.l1 = waverange[0]
    spider_params.l2 = waverange[1]

    p = pickle.load(open(path+"bestfit_"+model+".pic", "rb")) 
    d, m, par = p[0], p[1], p[2]            #data,  model, and best fit 

    if model == "zhang":
	spider_params.xi= par[d.par_order['xi']*d.nvisit]
        spider_params.T_n= par[d.par_order['T_n']*d.nvisit]
        spider_params.delta_T= par[d.par_order['T_n']*d.nvisit]

    elif model == "spherical":
        spider_params.degree=2
        spider_params.la0 = 0
        spider_params.lo0 = 0
  	sph0 = par[d.par_order['sph0']*d.nvisit]
  	sph1 = par[d.par_order['sph1']*d.nvisit]
  	sph2 = par[d.par_order['sph2']*d.nvisit]
  	sph3 = par[d.par_order['sph3']*d.nvisit]
        spider_params.sph = [sph0, sph1, sph2, sph3]
	print spider_params.sph

    elif model == "hotspot_t":
	spider_params.la0 = 0.
	spider_params.lo0 = 0.
	spider_params.size = 90.
	spider_params.spot_T = par[d.par_order['spot_T']*d.nvisit]
	spider_params.p_T = par[d.par_order['p_T']*d.nvisit]
        
    else:
        print "model not supported"

    t= spider_params.t0 + np.linspace(0, + spider_params.per,100)
    lc = spider_params.lightcurve(t)

    nla, nlo = 100, 100 
    #nla, nlo = 20, 20 

    las = np.linspace(-np.pi/2,np.pi/2,nla)
    los = np.linspace(-np.pi,np.pi,nlo)

    temp_map = True
    fluxes = []
    for la in las:
            line = []
            for lo in los:
                    flux = sp.call_map_model(spider_params,la,lo)
                    if temp_map == False:
                            line += [flux[0]]
                    else:
                            line += [flux[1]]
            fluxes += [line]
    
    fluxes = np.array(fluxes)
    #if model == "spherical": fluxes = convert_to_T(fluxes*np.pi*np.sqrt(2))

    vmax = 3900.
    vmin = 0.

    print "model, Tmin, Tmax = ", model, fluxes.min(), fluxes.max()
	

    im = plt.imshow(fluxes, aspect='auto', interpolation='bilinear', cmap=plt.cm.afmhot, alpha=0.8, vmax=vmax, vmin=vmin, extent = [-180,180, -90, 90])

    plt.gca().text(-160, 60, labels[i], bbox={'facecolor':'white', 'alpha':0.9, 'pad':5}, fontsize=10) 
    plt.xticks([-180, -90,  0, 90, 180])

    if col == 2: plt.ylabel("Latitude (deg)", labelpad =2)	
    if col == 3: plt.gca().set_yticklabels([])
    if row == 1: plt.xlabel("Longitude (deg)")
    if row == 0: plt.gca().set_xticklabels([])


ax = plt.subplot(gs[0:2, 4])
plt.gca().set_visible(False)
cb = fig.colorbar(im, alpha = 0.8, aspect = 30, fraction = 1) 
cb.set_label("Temperature (Kelvin)")

###############################
# plot best fits	      #
###############################

models = ["spherical", "hotspot_t", 'zhang']
path = "WFC3_best_fits/"
waverange = [1.1e-6, 1.7e-6]                    #WFC3
colors = ['#02ccfe', 'blue', '#ff7855'] 
grays = ['0.7', '0.5', '0.6']
labels = ["Sph. harmonics",  "Two temp.", "Kinematic"] 
#linestyle = ['-', '-.', ':']
linestyle = ['-', '-', '-']

phasebins = np.linspace(0., 1., 50)
nbins = len(phasebins) - 1

ax = plt.subplot(gs[0:2, 0])

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

	plt.errorbar(bin_average, (bindata-1.)*1e3*dilution, yerr = binsigma*1e3, linestyle='none', marker = '.', color = colors[i], ecolor = 'k', zorder = 10, alpha = 0.8)

        plt.plot(phase, (bestfit-1.)*1e3*dilution, color = colors[i], label = labels[i], alpha = 1.0, linestyle = linestyle[i])

plt.title("WFC3 broadband phase curve")
plt.legend(loc = 'upper left', frameon=False, fontsize=11)
plt.xlabel("Orbital phase")
plt.ylabel("Planet-to-star flux (ppt)")
plt.ylim(-0.2, 2)
plt.xlim(0,1)

plt.savefig('hst_model_comparison.pdf')
