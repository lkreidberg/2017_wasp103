import numpy as np
import scipy.interpolate as spi
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
import spiderman as sp
import pickle
from lc_fit import Model, LightCurveData

mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.size'] = 14
mpl.rcParams['lines.linewidth'] = 1.4

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
			if relative > 0: T[i,j] = get_T(relative)
			else: 
                                print relative
                                T[i,j] = 0.
	return T

def initialize_params(spider_params):
    spider_params.n_layers= 5
    spider_params.t0= 200               # Central time of PRIMARY transit [days]
    spider_params.per= 0.81347753       # Period [days]
    spider_params.a_abs= 0.01526        # The absolute value of the semi-major axis [AU]
    spider_params.inc= 82.33            # Inclination [degrees]
    spider_params.ecc= 0.0              # Eccentricity
    spider_params.w= 90                 # Argument of periastron
    spider_params.rp= 0.111            # Planet to star radius ratio
    spider_params.a= 4.855              # Semi-major axis scaled by stellar radius
    spider_params.p_u1= 0               # Planetary limb darkening parameter
    spider_params.p_u2= 0               # Planetary limb darkening parameter
    spider_params.T_s = 6110.		# Temperature of the star
    return spider_params

models = ["spherical", "zhang", "hotspot_t"]
labels = ["Spherical harmonics", "Kinematic", "Two temperature"] 
path = "WFC3_best_fits/"
waverange = [1.1e-6, 1.7e-6]                    #Spitzer Ch 2

fig =plt.figure(figsize = (6,6))

for i, model in enumerate(models):
    spider_params = sp.ModelParams(brightness_model=model)
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
    #print model
    #spider_params.square_plot()
    #plt.show()

    nla, nlo = 100, 100 
    #nla, nlo = 20, 20 

    las = np.linspace(-np.pi/2,np.pi/2,nla)
    los = np.linspace(-np.pi,np.pi,nlo)

    temp_map = True
    fluxes = []
    for la in las:
            row = []
            for lo in los:
                    flux = sp.call_map_model(spider_params,la,lo)
                    if temp_map == False:
                            row += [flux[0]]
                    else:
                            row += [flux[1]]
            fluxes += [row]
    
    fluxes = np.array(fluxes)
    if model == "spherical": fluxes = convert_to_T(fluxes*np.pi*np.sqrt(2))

    vmax = 4000.
    vmin = 0.

    print "model, Tmin, Tmax = ", model, fluxes.min(), fluxes.max()
	
    plt.subplot(311+i)

    #plt.imshow(fluxes, aspect='auto', interpolation='bilinear', cmap=plt.cm.YlOrRd_r, alpha=0.9, vmax=vmax, vmin=vmin, extent = [-180,180, -90, 90])
    im = plt.imshow(fluxes, aspect='auto', interpolation='bilinear', cmap=plt.cm.plasma, alpha=0.9, vmax=vmax, vmin=vmin, extent = [-180,180, -90, 90])

    plt.gca().text(-170, 50, labels[i], bbox={'facecolor':'white', 'alpha':0.9, 'pad':5}, fontsize=12)

    if i == 2: 
	plt.xticks([-180, -120,  -60,  0, 60, 120, 180])
	plt.xlabel("Longitude (degrees)")
    else: plt.xticks([])
    if i == 1: plt.ylabel("Latitude (degrees)")	


fig.subplots_adjust(right=0.75)
cbar_ax = fig.add_axes([0.8, 0.25, 0.05, 0.5])
cb = fig.colorbar(im, cax=cbar_ax)
cb.set_label("Temperature (Kelvin)")

plt.savefig('hst_maptemp.pdf')

