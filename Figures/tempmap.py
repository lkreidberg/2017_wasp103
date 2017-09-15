import numpy as np
import scipy.interpolate as spi
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
import spiderman as sp

mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.size'] = 14
mpl.rcParams['lines.linewidth'] = 1.4

#constants (all SI units)
h = 6.626e-34		#Planck constant
k = 1.38e-23		#Boltzmann constant
c = 3.0e8		#speed of light

#return blackbody spectral radiance as a function of wavelength lambda (in m) and temperature T (Kelvin)
def bb(l, T):
	return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T)-1)))

def get_T(relative):
	l = 6.e-6
	Ts = 4520.
	T = np.arange(300, 1600, 10)
	#rel = bb(l, T)/bb(l, Ts)*0.1594**2
	rel = bb(l, T)/bb(l, Ts)*0.06**2
	diff = np.abs(rel - relative)
	return T[diff == np.min(diff)]

def convert_to_T(flux):
	Ts = 4520
	l = 6.e-6
	n, m = flux.shape
	T = np.zeros_like(flux)
	for i in range(m):
		for j in range(n):
			relative = flux[i,j]
			T[i,j] = get_T(relative)
	return T

GCM = False

if GCM ==False:
	spider_params = sp.ModelParams(brightness_model="spherical")
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
#	spider_params.xi= 0.3       # Ratio of radiative to advective timescale
	#spider_params.T_n= 700	# Temperature of nightside
	#spider_params.delta_T= 900 # Day-night temperature contrast
	#spider_params.T_s = 4500    # Temperature of the star

	spider_params.degree=2
	spider_params.la0 = 0
	spider_params.lo0 = 0
	spider_params.sph = [2.75677373e-03, 2.33711161e-04, 1.88579054e-23, 1.21500964e-03]

	spider_params.l1 = 5.e-6       # The starting wavelength in meters
	spider_params.l2 = 12.e-6       # The ending wavelength in meters

	t= spider_params.t0 + np.linspace(0, + spider_params.per,100)
	lc = spider_params.lightcurve(t)


	nla, nlo = 100, 100

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
	fluxes = convert_to_T(fluxes)


else:
	d = np.genfromtxt("PTprofiles-WASP-43b-CLEAR-Solar-comp-C32.dat", delimiter = ',')
	pres = 0.02083
	ind = d[:,3] == pres
	fluxes = d[ind, 4]
	fluxes = fluxes.reshape(64,32).T



###############################################################################
# wavelength vs phase diagram
###############################################################################

tbmap9 = fluxes

vmax = tbmap9.max()
vmin = tbmap9.min()
plt.figure(5, figsize=(10,5))
plt.clf()
a = plt.axes([0.00,0.0,1.0,1.0])
plt.imshow(tbmap9 , aspect='auto', interpolation='bilinear', cmap=plt.cm.YlOrRd_r, alpha=0.9, vmax=vmax, vmin=vmin)
plt.savefig('map.png')

######################################
# Combined plot
######################################
#Tick marks
newxticks = np.linspace(-180,135, 8).astype(int)
oldxticks = np.linspace(0, 28, 8)
newyticks = np.arange(1.0, 1.71, 0.1)
oldyticks = np.linspace(0, 37, len(newyticks))

plt.figure(figsize=(4.5,3))
plt.clf()
vmin = fluxes.min()
vmax = fluxes.max()
ymin = -0.005
ymax =  0.080

#Temperature
m = Basemap(projection='robin',lon_0=0,resolution='c',celestial=False)
m.warpimage(image='map.png')
m.drawparallels(np.arange(-90.,120.,30.),labels=[False,True,False,False],labelstyle="+/-",size=12)
m.drawmeridians(np.arange(0.,360.,60.),labels=[False,False,False,True],labelstyle="+/-",size=12, rotation=30)
plt.savefig('tmap.pdf')

plt.clf()
fig = plt.figure(figsize=(4.5, 3))
ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.08])
cmap = mpl.cm.YlOrRd_r
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap,
                                norm=norm,
                                orientation='horizontal')
cb1.set_label('Temperature (Kelvin)', fontsize=18)


#COLORBAR
#for t in cb.ax.get_yticklabels():
#     t.set_fontsize(14)

plt.savefig('colorbar.pdf')
