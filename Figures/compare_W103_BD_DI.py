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
star_flux = star_flux/1.e3/star_wave	#[W/m^2/micron]
star_flux = star_flux*2.1074e20		#scale by r^2/d^2, w/ rstar = 1.436 Rs, d = 470 pc (from exoplanet.eu)

depth = 0.01195				#transit depth
scale = 1.e5

labels = ["nightside", "quadrature", "dayside"]

for i, f in enumerate(files):
	ax = plt.subplot(gs[i,0])

	w = np.linspace(1.0, 1.8, 100)*1.e-6
	T = Ts[i]
	flux = blackbody(w, T)/scale
	plt.plot(w*1e6, flux)
    	plt.plot(w*1e6, flux, color = 'w', linewidth=0., label = labels[i])

	d = np.genfromtxt(path+files[i]) 
	w = d[:,0]
	fs = np.interp(w, star_wave, star_flux)
	fp = d[:,1]*fs/depth/scale		#divide by depth to account for difference in rp vs rs
	fp_err = d[:,2]*fs/depth/scale
	
	plt.errorbar(w, fp, yerr = fp_err, fmt = '.k', zorder=100)
	
    	#plt.plot(w, fp, linewidth = 0., label = str(Ts[i])+ " K")
	plt.legend(loc = 'upper right', handlelength = 0., frameon=False, fontsize=9)

	plt.ylabel("W/m$^2$/$\mu$m")
	plt.xlim(1.0,1.8)
	plt.ylim(ylo[i], yhi[i])

	if i != 3: plt.gca().set_xticks([])
	if i == 3: plt.xlabel("Wavelength (microns)")
	if i == 1: plt.gca().set_yticks([5, 10, 15])
	if i == 0 : plt.title("WASP-103b")	


#Brown dwarf spectra
"""temp object
1918 +- 72  2MASSJ13204427p0409045_spec_app.txt
1988 +- 75  2MASSIJ2104149m103736_spec_app.txt
1995 +- 77  2MASSWJ1506544p132106_spec_app.txt
2212 +- 55  2MASPJ0345432p254023_spec_app.txt
2275 +- 56  LHS2924_spec_app.txt
2585 +- 46  2MASSJ03205965p1854233_spec_app.txt
2658 +- 121 LHS2021_spec_app.txt
2817 +- 59  2MASSJ03510004m0052452_spec_app.txt
2873 +- 75  2MASSJ00034227m2822410_spec_app.txt"""

path = "BD_spectra_original/"
files = ["2MASSJ13204427p0409045_spec_app.txt", "LHS2924_spec_app.txt", "2MASSJ00034227m2822410_spec_app.txt"]
names = ["1320+0409", "1428+3310", "0003-2822"]

Ts = [1918., 2275., 2873.]
distance = np.array([30.96, 11.01, 38.91])			#parsecs
radius = np.array([1.01, 1.06, 1.32])				#jupiter radii
unit_conversion = (distance/radius)**2*1.863e17*10.		#(r/d)^2 in dimensionless units *10 to get in W/m^2/micron from ergs/A/cm^2/s

for i, f in enumerate(files):

    ax = plt.subplot(gs[i,1])

    w = np.linspace(1.0, 1.8, 100)*1.e-6
    T = Ts[i]
    flux = blackbody(w, T)/scale
    plt.plot(w*1e6, flux)
    plt.plot(w*1e6, flux, color = 'w', linewidth=0., label = names[i])

    d = np.genfromtxt(path + f)
    plt.plot(d[:,0], d[:,1]*unit_conversion[i]/scale)

    #plt.plot(d[:,0], d[:,1], linewidth = 0., label = str(Ts[i])+ " K")
    plt.legend(loc = 'upper right', handlelength = 0., frameon=False, fontsize=9)
#    plt.gca().text(1, 1, str(Ts[i]))	#doesn't work bc it goes on last frame

    #plt.ylabel("W/m$^2$/$\mu$m")

    plt.xlim(1.0,1.8)
    plt.ylim(ylo[i], yhi[i])

    if i != 3: plt.gca().set_xticks([])
    if i == 3: plt.xlabel("Wavelength (microns)")
    if i == 1: plt.gca().set_yticks([5, 10, 15])
    if i == 0 : plt.title("Brown dwarfs")	


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
    plt.plot(d[:,0], d[:,1]*1.863e17*1.e2/guess)
    plt.plot(d[:,0], d[:,1]*1.863e17*1.e2/guess, 'w', linewidth=0., label = name[i])

    plt.legend(loc = 'upper right', handlelength = 0., frameon=False, fontsize=9)

    if i == 2: plt.ylim(ylo[i], yhi[i])
    plt.xlim(1.0,1.8)
    #plt.gca().set_yticks([])
    if i != 3: plt.gca().set_xticks([])
    if i == 3: plt.xlabel("Wavelength (microns)")
    if i == 0: plt.title("Imaged planets")

#plt.tight_layout()
plt.savefig("compare_W103_BD_DI.pdf")
#plt.show()
