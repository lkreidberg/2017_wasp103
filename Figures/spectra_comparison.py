import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy import modeling
import seaborn as sns

sns.set_context("paper")
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})

def get_blackbody_model(model_w, data_w, data, T):
    unit_conversion = 1.2533e19                                        #(10 pc/jupiter radius)^2

    model = modeling.blackbody.blackbody_lambda(model_w*1.e4, T)       #wavelength in angstroms
    model = model*4.*np.pi*1.e4*model_w/unit_conversion/np.pi          #multiplies by 4 pi steradians * 10000. angstroms/micron * micron
    #see chapter 1 of Rybicki & Lightman about the factor of pi to convert bb intensity to flux!
    scale = rescale(data_w, data, model_w, model)                                #rescales model to account for uncertainty in radius of the object
    return  model*scale
    

def rescale(data_w, data, model_w, model):
    model = np.interp(data_w, model_w, model)
    ind = (data > 1.e-20)
    data_w, data, model = data_w[ind], data[ind], model[ind]
    scale = np.linspace(0.05, 15, 1000.)
    chi2best, scalebest = 1.e10, -1.
    for s in scale:
        chi2 = np.sum((data - model*s)**2)
        if chi2 < chi2best:
                chi2best = chi2
                scalebest = s
    return scalebest 

gs = gridspec.GridSpec(3, 3, hspace=0.1, wspace=0.3)

#plotting parameters
C = 1.e10                                               #normalization constant for plotting
model_w = np.linspace(1, 1.75, 20)

#WASP-103b spectra
path = "W103_spectra/"
Ts = [1922, 2408, 2982]
files = ["phase_0.10_wfc3.txt", "phase_0.25_wfc3.txt", "phase_0.5_wfc3.txt"]

star = np.genfromtxt(path+"wasp103_sed_fluxes.out")
star_wave = star[:,0]	                        #[microns]
star_flux = star[:,1]	                        #[ergs/s/cm^2]	
star_flux = star_flux*(470/10.)**2		#convert to absolute flux (10 pc)

labels = ["nightside", "quadrature", "dayside"]
ylo = [0., 0.9, 2.8]
yhi = [1.15, 2.3, 4.3]

plt.figure(figsize = (7.5, 3.5))
#data_w, data, data_err; C; model_w, model

for i, f in enumerate(files):
    ax = plt.subplot(gs[i,0])

    #plot data
    d = np.genfromtxt(path+files[i]) 
    stellar_spectrum = np.interp(d[:,0], star_wave, star_flux)
    data_w = d[:,0]
    data = d[:,1]*stellar_spectrum
    data_err = d[:,2]*stellar_spectrum
    
    plt.errorbar(data_w, data*C, yerr = data_err*C, fmt = '.k', zorder=100)
    plt.plot(data_w, data*C, linewidth = 0., label = str(Ts[i])+ " K")
    
    #fit blackbody model
    model = get_blackbody_model(model_w, data_w, data, Ts[i])
    plt.plot(model_w, model*C, color = '0.5', linestyle = 'dotted', zorder  = -1)

    #labels and axes
    plt.legend(loc = 'upper right', handlelength = 0., frameon=False)

    if i == 1: plt.ylabel("erg/s/cm$^2$ ($\\times10^{-10}$)")
    plt.xlim(1.0,1.75)
    plt.ylim(ylo[i], yhi[i])

    if i != 2: plt.gca().set_xticklabels([])
    if i == 0 : plt.title("WASP-103b")	

#Brown dwarf spectra

path = "BD_spectra/"
files = ["For_Laura1320+0409 (L3:) SED.txt", "0428-2253 (L0.5) SED.txt", "For_Laura0003-2822 (M7.5) SED.txt"]
names = ["1320+0409", "0024-0158", "0003-2822"]
Ts = [1875, 2429, 2889]
distance = np.array([30.96, 25.99,38.91])
radius = np.array([1.01, 1.1,  1.33])
ylo = [0.05, 0.6, 2.0]
yhi = [0.5, 1.45, 4.1]

for i, f in enumerate(files):
    ax = plt.subplot(gs[i,1])

    d = np.genfromtxt(path + f)

    ind = (d[:,0]>1.1)&(d[:,0]<1.75)
    d = d[ind]

    data_w = d[:,0]
    data =  d[:,1]*1e4*d[:,0]
    plt.plot(data_w, data*C, color = 'k')             #original flux in ergs/cm^2/s/Angstrom (so multiply by wavelength in angstroms)
    plt.plot(data_w, data*C, linewidth = 0., label = str(Ts[i])+ " K")

    #fit blackbody model
    model = get_blackbody_model(model_w, data_w, data, Ts[i])
    plt.plot(model_w, model*C, color = '0.5', linestyle = 'dotted', zorder  = -1)

    #labels and axes
    plt.legend(loc = 'upper right', handlelength = 0., frameon=False)
    plt.xlim(1.0,1.75)
    plt.ylim(ylo[i], yhi[i])

    if i != 2: plt.gca().set_xticklabels([])
    if i == 2: plt.xlabel("Wavelength (microns)")
    if i == 0 : plt.title("Brown dwarfs")	

#Directly imaged spectra

path = "Direct_Image_Spectra/"
files = ["1RXSJ160929.1-210524b_Kreidberg.txt", "J1610-1913B_Kreidberg.txt", "TWA22A_Kreidberg.txt"]
name = ["1RXS", "J1610", "TWA22A"]
Ts = [1800, 2400, 3000] 
ylo = [0., 11., 11.]
yhi = [1.3, 23., 22.]

for i, f in enumerate(files):
    ax = plt.subplot(gs[i,2])

    d = np.genfromtxt(path + f)
    ind = (d[:,0]>1.1)&(d[:,0]<1.7)
    d = d[ind]
    data = d[:,1]*1.e3*d[:,0]
    data_w = d[:,0]
    
    plt.plot(data_w, data*C, color = 'k')

    #fit blackbody model
    model = get_blackbody_model(model_w, data_w, data, Ts[i])
    plt.plot(model_w, model*C, color = '0.5', linestyle = 'dotted', zorder  = -1)
    plt.plot(data_w, data*C, linewidth = 0., label = str(Ts[i])+ " K")

    #labels and axes
    plt.legend(loc = 'upper right', handlelength = 0., frameon=False)

    plt.xlim(1.0,1.75)
    plt.ylim(ylo[i], yhi[i])

    if i != 2: plt.gca().set_xticklabels([])
    if i == 0: plt.title("Imaged planets")

#plt.tight_layout()
plt.savefig("spectra_comparison.pdf")
#plt.show()
