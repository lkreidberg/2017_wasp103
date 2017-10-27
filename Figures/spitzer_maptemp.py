import numpy as np
import scipy.interpolate as spi
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
import spiderman as sp

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
	l = 4.5e-6
	Ts = 6110.
	rprs = 0.11092
	T = np.arange(10, 4000, 10)
	rel = bb(l, T)/bb(l, Ts)*rprs**2
	diff = np.abs(rel - relative)
	return T[diff == np.min(diff)]

def convert_to_T(flux):
	Ts = 6110. 
	l = 4.5e-6
	n, m = flux.shape
	T = np.zeros_like(flux)
	for i in range(m):
		for j in range(n):
			relative = flux[i,j]
			if relative > 0: T[i,j] = get_T(relative)
			else: T[i,j] = 0.
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

models = ['zhang', 'hotspot_t', "spherical"] 
labels = ["Kinematic", "Two Temperature", "Spherical harmonics"] 

plt.figure(figsize = (6,6))
channel = 2

for i, model in enumerate(models):
    spider_params = sp.ModelParams(brightness_model=model, thermal = True)
    spider_params = initialize_params(spider_params)
    if channel == 2: 
	#spider_params.filter = '/Users/lkreidberg/Desktop/Util/Throughput/spitzer_irac_ch2_electrons_per_energy.txt'
	spider_params.filter = '/Users/lkreidberg/Desktop/Util/Throughput/spitzer_irac_ch2.txt'
	spider_params.l1 = 3.9e-6
	spider_params.l2 = 5.1e-6
    elif channel == 1: 
	#spider_params.filter = '/Users/lkreidberg/Desktop/Util/Throughput/spitzer_irac_ch1_electrons_per_energy.txt'
	spider_params.filter = '/Users/lkreidberg/Desktop/Util/Throughput/spitzer_irac_ch1.txt'
	spider_params.l1 = 3.0e-6
	spider_params.l2 = 4.1e-6

    #print "still using filters that assume SPIDERMAN doesn't correct to units of energy -- update when you update spiderman"

    if model == "zhang":
	if channel == 2:
		spider_params.xi= 1.2977e-1       	# Ratio of radiative to advective timescale
		spider_params.T_n= 1703.		# Temperature of nightside
		spider_params.delta_T= 2606.		# Day-night temperature contrast
	elif channel == 1:
		spider_params.xi= 3.6969e-1       	# Ratio of radiative to advective timescale
		spider_params.T_n= 1955.		# Temperature of nightside
		spider_params.delta_T= 1857.		# Day-night temperature contrast

    elif model == "spherical":
	if channel == 2:
		spider_params.degree=2
		spider_params.la0 = 0
		spider_params.lo0 = 0
		spider_params.sph = [7.354e+3, 1.4627e+2, 0.0, 2.618e+3]
	elif channel == 1:
		spider_params.degree=2
		spider_params.la0 = 0
		spider_params.lo0 = 0
		spider_params.sph = [7.427e+3, 2.6104e+2, 0.0, 1.906e+3]

    elif model == "hotspot_t":
	if channel == 2:
		spider_params.la0 = 0.
		spider_params.lo0 = 0.
		spider_params.size = 90.
		spider_params.spot_T = 3311.
		spider_params.p_T = 1361.
	elif channel ==1:
		spider_params.la0 = 0.
		spider_params.lo0 = 0.
		spider_params.size = 90.
		spider_params.spot_T = 3041.
		spider_params.p_T = 1429.
        
    else:
        print "model not supported"

    #t= spider_params.t0 + np.linspace(0, + spider_params.per,100)
    #lc = spider_params.lightcurve(t)

    #nla, nlo = 100, 100 
    nla, nlo = 20, 20 

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
    #if model == "spherical": fluxes = convert_to_T(fluxes*np.pi)

    vmax = 4000. 
    vmin = 1000.

    print "model, Tmin, Tmax = ", model, fluxes.min(), fluxes.max()
	
    plt.subplot(311+i)

    #plt.imshow(fluxes, aspect='auto', interpolation='bilinear', cmap=plt.cm.YlOrRd_r, alpha=0.9, vmax=vmax, vmin=vmin, extent = [-180,180, -90, 90])
    plt.imshow(fluxes, aspect='auto', interpolation='bilinear', cmap=plt.cm.plasma, alpha=0.9, vmax=vmax, vmin=vmin, extent = [-180,180, -90, 90])

    plt.gca().text(-170, 50, labels[i], bbox={'facecolor':'white', 'alpha':0.9, 'pad':5}, fontsize=12)

    if i == 2: 
	plt.xticks([-180, -120,  -60,  0, 60, 120, 180])
	plt.xlabel("Longitude (degrees)")
    else: plt.xticks([])
    if i == 1: plt.ylabel("Latitude (degrees)")	

    #if i == 1: plt.colorbar()

plt.savefig('spitzer_maptemp.pdf')

