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

def proj_area(phi, inc):
        R0 = 1.07e8             #planet radius in m
        p = 0.00117             # planet-to-star mass ratio
        r = 2.9695e9            #orbital separation in m
        kn = 0.653              #from table B5 of Leconte et al. 2011
        n = 1.005
        qn = kn*(1.- n/5.)      #n 


        alpha1 = 2.5*qn*R0**3/(p*r**3)
        alpha2 = -1.25*qn*R0**3/(p*r**3)
        alpha3 = alpha2

        a1 = R0*(1+alpha1)
        a2 = R0*(1+alpha2)
        a3 = R0*(1+alpha3)

        return  np.sqrt(a3**2*np.sin(inc)**2*(a1**2*np.sin(phi)**2+a2**2*np.cos(phi)**2)+ a1**2*a2**2*np.cos(inc)**2)/(a2*a3)


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

models = ["hotspot_t", "zhang", "spherical"] 
labels = ["Two Temp.", "Kinematic", "Sph. Harmonics"] 
path = "WFC3_best_fits/"
waverange = [1.1e-6, 1.7e-6]                    #Spitzer Ch 2

fig =plt.figure(figsize = (10,3.33))
gs = gridspec.GridSpec(3, 5, width_ratios = [3, 0.1, 1,1,0.1], height_ratios = [1,1, 0.01]) 
ncol = 2

#cmap = plt.cm.afmhot
#cmap = plt.cm.plasma
#cmap = plt.cm.magma
#cmap = plt.cm.viridis
#cmap = plt.cm.Purples_r
cmap = plt.cm.inferno

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
        print "sph", spider_params.sph
	#print spider_params.sph

    elif model == "hotspot_t":
	spider_params.la0 = 0.
	spider_params.lo0 = 0.
	spider_params.size = 90.
	spider_params.spot_T = par[d.par_order['spot_T']*d.nvisit]
	spider_params.p_T = par[d.par_order['p_T']*d.nvisit]

    elif model == "sincos":
	print "AHHH"
        
    else:
        print "model not supported"

    t= spider_params.t0 + np.linspace(0, + spider_params.per,100)
    lc = spider_params.lightcurve(t)

    nla, nlo = 100, 100 
    #nla, nlo = 10, 10 

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
#			    print la, lo, flux[1]
            fluxes += [line]
    
    fluxes = np.array(fluxes)
    #if model == "spherical": fluxes = fluxes/np.pi 

    vmax = 3900.
    vmin = 0.

    print "model, Tmin, Tmax = ", model, fluxes.min(), fluxes.max()
	

    im = plt.imshow(fluxes, aspect='auto', interpolation='bilinear', cmap=cmap, alpha=0.9, vmax=vmax, vmin=vmin, extent = [-180,180, -90, 90])

    plt.gca().text(-160, 60, labels[i], bbox={'facecolor':'white', 'alpha':0.9, 'pad':5}, fontsize=10) 
    plt.xticks([-180, -90,  0, 90, 180])

    if col == 2: plt.ylabel("Latitude (deg)", labelpad =2)	
    if col == 3: plt.gca().set_yticklabels([])
    if row == 1: plt.xlabel("Longitude (deg)")
    if row == 0: plt.gca().set_xticklabels([])

#plots GCM
GCM = np.genfromtxt("GCM_From_Vivien/PTprofiles-WASP-103b-TiO-fix-3-Drag3-NEW-OPA-nit-1036800.dat", delimiter = ",")
ax = plt.subplot(gs[1,3])
Pref = 0.11542					#reference pressure in bars
#Pref = 0.010131
ind = GCM[:,3] == Pref
lat, lon = GCM[ind,1], GCM[ind,2]
Ts = GCM[ind,4]

fluxes = np.reshape(Ts, (-1, 30))
print "max, min T", np.max(Ts), np.min(Ts)
im = plt.imshow(fluxes.T, aspect='auto', interpolation='bilinear', cmap=cmap, alpha=0.9, vmax=vmax, vmin=vmin, extent = [-180,180, -90, 90])
plt.gca().text(-160, 60, "GCM", bbox={'facecolor':'white', 'alpha':0.9, 'pad':5}, fontsize=10) 
plt.xticks([-180, -90,  0, 90, 180])
plt.gca().set_yticklabels([])
plt.xlabel("Longitude (deg)")


#formatting
ax = plt.subplot(gs[0:2, 4])
plt.gca().set_visible(False)
cb = fig.colorbar(im, alpha = 0.8, aspect = 30, fraction = 1) 
cb.set_label("Temperature (Kelvin)")


###############################
# plot best fits	      #
###############################

models = ["hotspot_t", "zhang", "spherical", "sincos"]
labels = ["Two temp.", "Kinematic", "Sph. harmonics", "Sinusoid"] 
path = "WFC3_best_fits/"
waverange = [1.1e-6, 1.7e-6]                    #WFC3
colors = ['#82cafc', '#ff7855', 'purple', 'red'] 
grays = ['0.7', '0.5', '0.6']
#linestyle = ['-', '-.', ':']
linestyle = ['-', '-', '-', '-']

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

	#plt.errorbar(bin_average, (bindata-1.)*1e3*dilution, yerr = binsigma*1e3, linestyle='none', marker = '.', color = colors[i], ecolor = 'k', zorder = 10, alpha = 0.8)
	plt.plot(bin_average, (bindata-1.)*1e3*dilution,  linestyle = 'none', marker ='.', color = colors[i], zorder = 10, alpha = 0.8)
        pickle.dump( [bin_average, (bindata-1.), binsigma], open(model+".p", "wb" ) )

        plt.plot(phase, (bestfit-1.)*1e3*dilution, color = colors[i], label = labels[i], alpha = 0.5, linestyle = linestyle[i])


#plot GCM
inc = 1.523
#gcms = ["./Vivien_models2/PC-TiO-NoClouds-Drag3-NEW-OPA-NEW-PT.dat"]
gcms = ["./Vivien_models2/PC-TiO-NoClouds-Drag1-NEW-OPA-NEW-PT.dat"]
#gcms = ["./Vivien_models2/PC-TiO-NoClouds-NEW-OPA-NEW-PT.dat"]
correction_factor = 1.099
#plot gcms
for i, g in enumerate(gcms):
    model = np.genfromtxt(g, delimiter = ',') 
    phi = model[:,0]/360.*2.*np.pi
    area = proj_area(phi, inc)
    plt.plot(model[:,0]/360. +0.5, model[:,5]*1e3*area*correction_factor, label = "GCM", color = '0.7', zorder = -20, linestyle = '--') 


plt.title("WFC3 broadband phase curve")
plt.legend(loc = 'upper left', frameon=False, fontsize=11)

ax = plt.gca()
handles,labels = ax.get_legend_handles_labels()

handles = [handles[2], handles[1], handles[0], handles[3], handles[4]]
labels = [labels[2], labels[1], labels[0], labels[3], "GCM"]

ax.legend(handles,labels,loc='upper left', frameon = False, fontsize = 11)


plt.xlabel("Orbital phase")
plt.ylabel("Planet-to-star flux ($\\times10^{-3}$)")
plt.ylim(-0.2, 2)
plt.xlim(0,1)

#gs.tight_layout(fig, rect = [0., 0., 1., 1.])
plt.savefig('hst_model_comparison.pdf')


#plt.show()
