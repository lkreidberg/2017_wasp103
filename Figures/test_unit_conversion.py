import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pyfits
from astropy import units as u

def blackbody(l,T): 
	h = 6.62607e-34   #J/s
	c = 2.997925e8       #m/s
	k = 1.38065e-23     #J/K
	#return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T))-1))			#[W/sr/m^3]
	return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T))-1))*np.pi/1.e6	#[W/m^2/micron]

gs = gridspec.GridSpec(3, 3, hspace=0.1, wspace=0.3)

ylo = [0, 5, 12]
yhi = [6, 15, 33]

#WASP-103b spectra
path = "W103_spectra/"
Ts = [1911., 2408., 2921.] 
files = ["phase_0.10_wfc3.txt", "phase_0.25_wfc3.txt", "phase_0.5_wfc3.txt"]

star = np.genfromtxt(path+"wasp103_sed_fluxes.out")
star_wave = star[:,0]	#[microns]
star_flux = star[:,1]	#[ergs/s/cm^2]	
star_flux = star_flux*(470/10.)**2		#convert to absolute flux (10 pc)


labels = ["nightside", "quadrature", "dayside"]

for i, f in enumerate(files):
	ax = plt.subplot(gs[i,0])

	w = np.linspace(1.0, 1.8, 100)*1.e-6
	T = Ts[i]

	d = np.genfromtxt(path+files[i]) 
	w = d[:,0]
	fs = np.interp(w, star_wave, star_flux)
	fp = d[:,1]*fs
	fp_err = d[:,2]*fs
	
	plt.errorbar(w, fp, yerr = fp_err, fmt = '.k', zorder=100)
	
    	#plt.plot(w, fp, linewidth = 0., label = str(Ts[i])+ " K")
	plt.legend(loc = 'upper right', handlelength = 0., frameon=False, fontsize=9)

	plt.ylabel("W/m$^2$/$\mu$m")
	plt.xlim(1.0,1.8)
	#plt.ylim(ylo[i], yhi[i])

	if i != 3: plt.gca().set_xticks([])
	if i == 3: plt.xlabel("Wavelength (microns)")
	if i == 0 : plt.title("WASP-103b")	

#Directly imaged spectra
"""1RXSJ160929.1-210524b   L2gamma  1800K
ABPicB  L0gamma 1700K
CD-352722B      L3      1800K
GSC08047-00232B M9.5    2200K
HIP78530B       M9.5    2800K
2MASS J1610-1913B       M9      2400K
PZTelB  M7      2700K
TWA11C  M5      3100K
TWA22A  M5      3000K
TWA22B  M5.5 2900K
USCOCTIO108B    M9.5    2300K"""


path = "Direct_Image_Spectra/"
files = ["1RXSJ160929.1-210524b_Kreidberg.txt", "J1610-1913B_Kreidberg.txt", "TWA22A_Kreidberg.txt"]
name = ["1RXS", "J1610", "TWA22A"]

guess = 1.e6
for i, f in enumerate(files):
    ax = plt.subplot(gs[i,2])
    d = np.genfromtxt(path + f)
    #plt.plot(d[:,0], d[:,1]*1.863e17*1.e2/guess)
    plt.plot(d[:,0], d[:,1]*1.e3*d[:,0])

    plt.legend(loc = 'upper right', handlelength = 0., frameon=False, fontsize=9)

    #if i == 2: plt.ylim(ylo[i], yhi[i])
    plt.xlim(1.0,1.8)
    #plt.gca().set_yticks([])
    if i != 3: plt.gca().set_xticks([])
    if i == 3: plt.xlabel("Wavelength (microns)")
    if i == 0: plt.title("Imaged planets")

#plt.tight_layout()
plt.show()
