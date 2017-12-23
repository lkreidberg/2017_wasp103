import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy import modeling

def rescale(data, model):
    ind = data == 0
    data = np.ma.array(data, mask = ind)
    scale = np.linspace(0.01, 10, 100.)
    chi2best, scalebest = 1.e10, -1.
    for s in scale:
        chi2 = np.sum((data - model*s)**2)
        if chi2 < chi2best:
                chi2best = chi2
                scalebest = s
    return scalebest


gs = gridspec.GridSpec(3, 3, hspace=0.1, wspace=0.3)
C = 1.e10                               #normalization constant for plotting

#WASP-103b spectra
path = "W103_spectra/"
#Ts = [1911., 2408., 2921.] 
Ts = [1922., 2408., 2982.]
files = ["phase_0.10_wfc3.txt", "phase_0.25_wfc3.txt", "phase_0.5_wfc3.txt"]

star = np.genfromtxt(path+"wasp103_sed_fluxes.out")
star_wave = star[:,0]	                        #[microns]
star_flux = star[:,1]	                        #[ergs/s/cm^2]	
star_flux = star_flux*(470/10.)**2		#convert to absolute flux (10 pc)


labels = ["nightside", "quadrature", "dayside"]

for i, f in enumerate(files):
	ax = plt.subplot(gs[i,0])

	d = np.genfromtxt(path+files[i]) 
        fs = np.interp(d[:,0], star_wave, star_flux)
	fp = d[:,1]*fs
	fp_err = d[:,2]*fs
	
        plt.errorbar(d[:,0], fp*C, yerr = fp_err*C, fmt = '.k', zorder=100)
	
    	#plt.plot(w, fp, linewidth = 0., label = str(Ts[i])+ " K")

        w = np.linspace(1, 2, 20)
        bb = modeling.blackbody.blackbody_lambda(w*1.e4, Ts[i])       #wavelength in angstroms
        bb = bb*4.*np.pi*1.e4*w        #multiplies by 4 pi steradians * 10000. angstroms/micron * micron
        scale = 1.2533e19  #(10 pc/jupiter radius)^2
        blackbody =  bb/scale/np.pi

        model = np.interp(d[:,0], w, blackbody)
        plt.plot(w, blackbody*rescale(fp,model)*C, color = '0.5', linestyle = 'dotted')

	plt.legend(loc = 'upper right', handlelength = 0., frameon=False, fontsize=9)

        if i == 1: plt.ylabel("erg/s/cm$^2$ [$\\times10^{-10}$]")
	plt.xlim(1.0,1.8)
	#plt.ylim(0.5*np.min(blackbody.value), 1.7*np.max(blackbody.value))

	if i != 2: plt.gca().set_xticklabels([])
	if i == 0 : plt.title("WASP-103b")	

#Brown dwarf spectra

path = "BD_spectra/"
files = ["For_Laura1320+0409 (L3:) SED.txt", "0428-2253 (L0.5) SED.txt", "For_Laura0003-2822 (M7.5) SED.txt"]
names = ["1320+0409", "0024-0158", "0003-2822"]
Ts = [1875., 2429., 2889.]
distance = np.array([30.96, 25.99,38.91])
radius = np.array([1.01, 1.1,  1.33])

normalization = radius**2*5.368e-20     #radius*(1 rjup/10 pc)**2        

for i, f in enumerate(files):

    ax = plt.subplot(gs[i,1])

    w = np.linspace(1, 2, 10)
    bb = modeling.blackbody.blackbody_lambda(w*1.e4, Ts[i])       #wavelength in angstroms
    bb = bb*4.*np.pi*1.e4*w        #multiplies by 4 pi steradians * 10000. angstroms/micron * micron
    plt.plot(w, bb*normalization[i]/np.pi*C, color = '0.5', linestyle = 'dotted')

    d = np.genfromtxt(path + f)
    plt.plot(d[:,0], d[:,1]*1e4*d[:,0]*C, color = 'k')             #original flux in ergs/cm^2/s/Angstrom (so multiply by wavelength in angstroms)

    #plt.plot(d[:,0], d[:,1], linewidth = 0., label = str(Ts[i])+ " K")
    plt.legend(loc = 'upper right', handlelength = 0., frameon=False, fontsize=9)
#    plt.gca().text(1, 1, str(Ts[i]))	#doesn't work bc it goes on last frame


    plt.xlim(1.0,1.8)
    #plt.ylim(ylo[i], yhi[i])

    if i != 2: plt.gca().set_xticks([])
    if i == 2: plt.xlabel("Wavelength (microns)")
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
Ts = [1800., 2400., 3000.] 

guess = 1.e6
for i, f in enumerate(files):
    ax = plt.subplot(gs[i,2])
    d = np.genfromtxt(path + f)
    ind = (d[:,0]>1.1)&(d[:,0]<1.7)
    d = d[ind]
    data = d[:,1]*1.e3*d[:,0]
    
    w = np.linspace(1, 2, 20)
    bb = modeling.blackbody.blackbody_lambda(w*1.e4, Ts[i])       #wavelength in angstroms
    bb = bb*4.*np.pi*1.e4*w                     #multiplies by 4 pi steradians * 10000. angstroms/micron * micron
    normalization = 1.2533e19*np.pi             #(10 pc/jupiter radius)^2
    bb = bb/normalization
    model = np.interp(d[:,0], w, bb)
    plt.plot(w, bb*rescale(data, model)*C, color = '0.5', linestyle = 'dotted') 

    plt.plot(d[:,0], data*C, color = 'k')
    plt.legend(loc = 'upper right', handlelength = 0., frameon=False, fontsize=9)

    plt.xlim(1.0,1.8)
    if i != 2: plt.gca().set_xticks([])
    if i == 0: plt.title("Imaged planets")

#plt.tight_layout()
plt.savefig("test_unit_conversion.pdf")
#plt.show()
