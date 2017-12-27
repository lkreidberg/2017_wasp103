import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy import modeling
import seaborn as sns

sns.set_context("talk")
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})

def rescale(data, model):
    ind = data > 1.e-20
    data, model = data[ind], model[ind]
    #data = np.ma.array(data, mask = ind)
    scale = np.linspace(0.01, 10, 1000.)
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

        w = np.linspace(1, 1.75, 20)
        bb = modeling.blackbody.blackbody_lambda(w*1.e4, Ts[i])       #wavelength in angstroms
        bb = bb*4.*np.pi*1.e4*w        #multiplies by 4 pi steradians * 10000. angstroms/micron * micron
        scale = 1.2533e19  #(10 pc/jupiter radius)^2
        blackbody =  bb/scale/np.pi

        model = np.interp(d[:,0], w, blackbody)
        blackbody = blackbody*rescale(fp,model)*C 
        plt.plot(w, blackbody, color = '0.5', linestyle = 'dotted', zorder  = -1)

	plt.legend(loc = 'upper right', handlelength = 0., frameon=False, fontsize=9)

        if i == 1: plt.ylabel("erg/s/cm$^2$ ($\\times10^{-10}$)")
	plt.xlim(1.0,1.75)
	plt.ylim(np.median(blackbody.value) - 0.6, np.median(blackbody.value) + 0.6)

	if i != 2: plt.gca().set_xticklabels([])
	if i == 0 : plt.title("WASP-103b")	

#Brown dwarf spectra

path = "BD_spectra/"
files = ["For_Laura1320+0409 (L3:) SED.txt", "0428-2253 (L0.5) SED.txt", "For_Laura0003-2822 (M7.5) SED.txt"]
names = ["1320+0409", "0024-0158", "0003-2822"]
Ts = [1875., 2429., 2889.]
distance = np.array([30.96, 25.99,38.91])
radius = np.array([1.01, 1.1,  1.33])

normalization = radius**2*5.368e-20/np.pi     #radius*(1 rjup/10 pc)**2        

for i, f in enumerate(files):

    ax = plt.subplot(gs[i,1])

    d = np.genfromtxt(path + f)

    ind = (d[:,0]>1.1)&(d[:,0]<1.8)
    d = d[ind]
    data =  d[:,1]*1e4*d[:,0]
    plt.plot(d[:,0], data*C, color = 'k')             #original flux in ergs/cm^2/s/Angstrom (so multiply by wavelength in angstroms)

    w = np.linspace(1, 1.75, 10)
    bb = modeling.blackbody.blackbody_lambda(w*1.e4, Ts[i])       #wavelength in angstroms
    bb = bb*4.*np.pi*1.e4*w        #multiplies by 4 pi steradians * 10000. angstroms/micron * micron
    bb = bb*normalization[i]
    model = np.interp(d[:,0], w, bb)
    blackbody = bb*rescale(data, model)*C 
    plt.plot(w, blackbody, color = '0.5', linestyle = 'dotted', zorder = -1)

    plt.legend(loc = 'upper right', handlelength = 0., frameon=False, fontsize=9)

    #plt.xlim(1.0,1.8)
    plt.xlim(1.0,1.75)
    plt.ylim(np.median(blackbody.value) - 0.6, np.median(blackbody.value) + 0.6)

    if i != 2: plt.gca().set_xticks([])
    if i == 2: plt.xlabel("Wavelength (microns)")
    if i == 0 : plt.title("Brown dwarfs")	

#Directly imaged spectra

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
    
    w = np.linspace(1, 1.75, 20)
    bb = modeling.blackbody.blackbody_lambda(w*1.e4, Ts[i])       #wavelength in angstroms
    bb = bb*4.*np.pi*1.e4*w                     #multiplies by 4 pi steradians * 10000. angstroms/micron * micron
    normalization = 1.2533e19*np.pi             #(10 pc/jupiter radius)^2
    bb = bb/normalization
    model = np.interp(d[:,0], w, bb)
    blackbody = bb*rescale(data, model)*C 
    plt.plot(w, blackbody, color = '0.5', linestyle = 'dotted', zorder = -1) 

    plt.plot(d[:,0], data*C, color = 'k')
    plt.legend(loc = 'upper right', handlelength = 0., frameon=False, fontsize=9)

    plt.xlim(1.0,1.75)
    plt.ylim(np.median(blackbody.value) - 0.6, np.median(blackbody.value) + 0.6)

    plt.xlim(1.0,1.75)
    if i != 2: plt.gca().set_xticks([])
    if i == 0: plt.title("Imaged planets")

#plt.tight_layout()
print "TODO:"
print "figure out what y axis scaling you want"
print "why is bb normalization off for middle temperature DI planet?"
print "STANDARDIZE notation"
plt.savefig("test_unit_conversion.pdf")
#plt.show()
