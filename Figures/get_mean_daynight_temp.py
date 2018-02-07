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


def initialize_params(spider_params):
    spider_params.n_layers= 5
    spider_params.t0= 200               # Central time of PRIMARY transit [days]
    spider_params.per= 0.9255       # Period [days]
    spider_params.a_abs= 0.01526        # The absolute value of the semi-major axis [AU]
    spider_params.inc= 87.3            # Inclination [deg]
    spider_params.ecc= 0.0              # Eccentricity
    spider_params.w= 90                 # Argument of periastron
    spider_params.rp= 0.1146            # Planet to star radius ratio
    spider_params.a= 2.999              # Semi-major axis scaled by stellar radius
    spider_params.p_u1= 0               # Planetary limb darkening parameter
    spider_params.p_u2= 0               # Planetary limb darkening parameter
    spider_params.T_s = 6110.		# Temperature of the star
    return spider_params

models = ["hotspot_t", "zhang", "spherical"]
labels = ["Two temp.", "Kinematic", "Sph. harmonics"] 

path = "WFC3_best_fits/"
waverange = [1.1e-6, 1.7e-6]        

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

    if model == "spherical": spider_params = sp.ModelParams(brightness_model=model, thermal=True, stellar_model = "w103_spectrum.txt")
    else: spider_params = sp.ModelParams(brightness_model=model, stellar_model = "w103_spectrum.txt")

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
#	print spider_params.sph

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

    nla, nlo = 200, 200 
    #nla, nlo = 300, 300 

    las = np.linspace(-np.pi/2,np.pi/2,nla)
    los = np.linspace(-np.pi,np.pi,nlo)

    dayside = []
    nightside = []
    Ts = np.zeros((nla, nlo))
    for i, la in enumerate(las):
            for j, lo in enumerate(los):
                    flux = sp.call_map_model(spider_params,la,lo)
                    Ts[i,j] = flux[1]
                    if np.abs(lo)<np.pi/2.: dayside.append(flux[1])
                    else: nightside.append(flux[1])

    
    dayside = np.array(dayside)
    nightside = np.array(nightside)
    
    Ts = np.array(Ts)
    #if model == "spherical": Ts, dayside, nightside = Ts/np.pi, dayside/np.pi, nightside/np.pi

    """lats, lons = np.meshgrid(las, los)
    ind = np.abs(lons)<np.pi/2."""

#    print "model, Tmin, Tmax = ", model, Ts.min(), Ts.max()
    print model, "daymu, nightmu, min, max", np.mean(dayside), np.mean(nightside), Ts.min(), Ts.max()
	

