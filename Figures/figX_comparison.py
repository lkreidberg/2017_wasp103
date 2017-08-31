import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pyfits
from astropy import units as u

def blackbody(l,T): 
	h = 6.626e-34   #J/s
	c = 3.0e8       #m/s
	k = 1.4e-23     #J/K
	return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T))-1))

gs = gridspec.GridSpec(4, 3, hspace=0.05, wspace=0.25)

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

path = "BD_spectra/"
files = ["2MASSJ13204427p0409045_spec_app.txt", "2MASPJ0345432p254023_spec_app.txt", "2MASSJ03205965p1854233_spec_app.txt", "2MASSJ00034227m2822410_spec_app.txt"]

w = np.linspace(0.3, 10, 100)*1.e-6
#T = 3130.
"""T = 2000.
flux = blackbody(w, T)
plt.plot(w*1e6, flux/1.e6*np.pi)
plt.gca().set_yscale('log')
plt.gca().set_xscale('log')
plt.show()"""

Ts = [1918., 2212., 2585., 2873.]

for i, f in enumerate(files):

    ax = plt.subplot(gs[i,0])

    w = np.linspace(0.5, 2, 100)*1.e-6
    T = Ts[i]
    flux = blackbody(w, T)
    plt.plot(w*1e6, flux/1.e6*np.pi)
    plt.plot(w*1e6, flux/1.e6*np.pi, color = 'w', linewidth=0., label = f[5:13])

    d = np.genfromtxt(path + f)
    plt.plot(d[:,0], d[:,1]*1.e21)

    plt.plot(d[:,0], d[:,1]*1.e21, linewidth = 0., label = str(Ts[i])+ " K")
    plt.legend(loc = 'lower right', handlelength = 0., frameon=False, fontsize=9)
#    plt.gca().text(1, 1, str(Ts[i]))	#doesn't work bc it goes on last frame

    ax.set_yscale('log')
    plt.ylabel("W/m$^2$/$\mu$m")

    plt.xlim(0.5,2)
    if i != 3: plt.gca().set_xticks([])
    if i == 3: plt.xlabel("Wavelength (microns)")
    if i == 0 : plt.title("Brown dwarfs")	


#Directly imaged spectra
"""2M1207 A     2500  100 3.5  0.5 2640
OTS 44          1700  100 3.5  0.5 2300
KPNO Tau 4      1700  100 3.5  0.5 2300
2M0141          1800  200 4.5  0.5 2200
TWA 5B          2500  100 4.0  0.5 2570
Cha 1109        1800  100 4.0  0.5 2300
Gl 417 B        1800  100 5.0  0.5 1660
AB Pic b        1800  200 4.5  0.5 2200
CT Cha b        2600  100 3.5  0.5 2700
DH Tau B        2400  100 3.5  0.5 2350
GSC8047 B       2200  100 4.0  0.5 2300
USco CTIO 108B  2300  100 4.0  0.5 2300
HR7329 B        2600  100 . . . 2570
2M0345          2400  100 5.0  0.5 2230
TWA22 A         3000  100 4.5  0.5 3125
TWA22 B         3000  100 4.5  0.5 3064"""

path = "DI_spectra/"
files = ["L0_R2000_ABPicb_J_0.025_norm_master_spectrum.fits", "M8.5_R2000_2M1207A_SINFONIspeclib_JHK.fits", "M6_R2000_TWA8B.fits", "M5_R1200_TWA11C.fits"]
name = ["L0", " M8.5", "M6", "M5"]

for i, f in enumerate(files):
    ax = plt.subplot(gs[i,1])
    p = pyfits.open(path + f)
    d = p[0].data 
    if i>1: d = d.T
    plt.plot(d[:,0], d[:,1])
    plt.plot(d[:,0], d[:,1], linewidth = 0., label = name[i])

    if i == 0: 
	p = pyfits.open(path + "L0_R2000_ABPicb_H+K_0.025_final_spectrum.fits")
	d = p[0].data 
        plt.plot(d[:,0], d[:,1])

    plt.legend(loc = 'lower left', handlelength = 0., frameon=False, fontsize=9)

    plt.xlim(0.5,2)
    plt.gca().set_yticks([])
    if i != 3: plt.gca().set_xticks([])
    if i == 3: plt.xlabel("Wavelength (microns)")
    if i == 0: plt.title("Imaged planets")

plt.savefig("comparison.pdf")
plt.show()
