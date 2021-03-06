import pickle
import numpy as np
import matplotlib.pyplot as plt
import os, glob
from lc_fit import Model, LightCurveData
#import read_spectra_from_Vivien as rs
import matplotlib.gridspec as gridspec
from pylab import *
from matplotlib import rc, ticker
import seaborn as sns
from astropy.convolution import Gaussian1DKernel, convolve
import scipy.stats as st

sns.set_context("notebook", font_scale = 1.5)
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def get_significance(chisq, dof):
        alpha = (1. - st.chi2.cdf(chisq, dof))/2.
        z = st.norm.ppf(1.-alpha)
        return z

def blackbody(l,T): 
	h = 6.62607e-34   #J/s
        c = 2.997925e8       #m/s
        k = 1.38065e-23     #J/K
        return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T))-1))                 #[W/sr/m^3]

def best_fit_bb(w, y, e, rprs):
    Ts = np.linspace(1000, 3600, 10000)
    chibest = 10000.
    Tbest = 0.  
    Tstar = 6110.
    resid = 0.

    w = np.array(w)

    w, y, e = w[0:10], y[0:10], e[0:10]     #WFC3 data only
    #w, y, e = w[10], y[10], e[10]           #spitzer ch1
    #w, y, e = w[11], y[11], e[11]           #spitzer ch2

    #get stellar spectrum
    g = Gaussian1DKernel(stddev=100)
    star = np.genfromtxt("W103-6110K-1.22Ms-1.22Rs.dat")       #W/m2/micron (column 1)
    star_bb = np.interp(w, star[:,0], convolve(star[:,1]*22423.,g))
    #star = np.genfromtxt("wasp103_sed_fluxes.out")
    #star_bb = np.interp(w, star[:,0], star[:,1])*1.e24/(w*np.pi*4.)
    chis = []

    for T in Ts:
        #model = blackbody(w*1.0e-6, T)/star_bb*rprs**2
        model = (np.pi/1.e6)*blackbody(w*1.0e-6, T)/star_bb*rprs**2
        chi2 = np.sum((y - model)**2/e**2)
        chis.append(chi2)
        if chi2 < chibest: chibest, Tbest, resid = chi2, T, (y - model)/e

    waves_hires = np.linspace(1.0, 5.0, 100)
    #star_bb_hires = np.interp(waves_hires, star[:,0], star[:,1])*1.e24/(waves_hires*np.pi*4.)
    star_bb_hires = np.interp(waves_hires, star[:,0], convolve(star[:,1]*22423., g))

    outfile = open("temperatures.txt", "a")

    idx = (np.abs(chis-(chibest+1.))).argmin()
    onesigma = Tbest - Ts[idx]
    #print '{0:d}'.format(int(Tbest)),  '{0:d}'.format(int(np.abs(onesigma))), '{0:0.1f}'.format(chibest/(len(y) - 1)),
    #print '{0:d}'.format(int(Tbest)),  '{0:d}'.format(int(np.abs(onesigma))), get_significance(chibest,(len(y) - 1))
    print '{0:d}'.format(int(Tbest)),  '{0:d}'.format(int(np.abs(onesigma))) 
    print>>outfile,  Tbest                                                                                                            
    outfile.close() 
    return waves_hires, np.pi/1.e6*blackbody(waves_hires*1.0e-6, Tbest)/star_bb_hires*rprs**2 

GCM = np.genfromtxt("Vivien_models2/SpectralPC-Phi-TiO-NoClouds-Drag1-NEW-OPA-NEW-PT.dat", delimiter = ",")

path  = "WFC3_best_fits/spec_fits/"		#zhang model
files = glob.glob(os.path.join(path, "bestfit*.pic"))		


phasebins = np.array([0.06, 0.15, 0.25, 0.35, 0.44])
phasebins2 = np.array([0.56, 0.65, 0.75, 0.85, 0.94]) 

phasebins = np.append(phasebins, 0.5)
phasebins = np.append(phasebins, phasebins2)

ymins = [0.04, 0.2, 0.5, 0.6,  0.6, 0.5, 0.2, 0.04]
ymaxs = [3, 4, 5, 6, 6, 5, 4, 3]

nspitz_chan = 2
spectra = np.zeros((len(files) + nspitz_chan, len(phasebins)-1, 2))
phases = np.zeros(len(phasebins) -1 )
waves, dilution = [], []
i = 0

for f in files:
	p = pickle.load(open(f, 'rb'))
	d, m = p[0], p[1]			#stores data and model into d & m
	
	waves.append(d.wavelength)
	dilution.append(d.dilution)

	ind = d.err < 9.0e7			#indices of outliers
	
	err = d.err[ind]/d.flux[ind]		#normalized per point uncertainty 
	phase = m.phase[ind]
	data_corr = m.data_corr[ind]

	for j in range(1, len(phasebins)):
		#calculates uncertainty for in-eclipse points
		ind1 = (phase>=0.45)&(phase<=0.55)
		sig1 = np.sqrt(np.sum(err[ind1]**2)/sum(ind1)**2)

		#calculates uncertainty for phase bin
		ind2 = (phase >= phasebins[j-1]) & (phase < phasebins[j])
		sig2 = np.sqrt(np.sum(err[ind2]**2)/sum(ind2)**2)
	#	print "std, se",  np.sqrt(np.sum(err[ind2]**2)/sum(ind2)**2), np.sqrt(np.sum(err[ind2]**2)/(sum(ind2)-1)**2)            #difference of 1 ppm for std vs ste

		spectra[i, j-1, 0] = np.mean(data_corr[ind2])
		spectra[i, j-1, 1] = np.sqrt(sig1**2 + sig2**2)		#quadrature add uncertainties from baseline and bin 
		phases[j-1] = (phasebins[j] + phasebins[j-1])/2.

	i += 1

#add Spitzer Ch 1
#f= "Ch1_best_fits/2017-10-11_20:25-zhang/bestfit.pic"
f= "Ch1_best_fits/2018-02-07_14:28-spherical/bestfit.pic"
waves.append(3.6)
dilution.append(0.1712)

p = pickle.load(open(f, 'rb'))
data_corr = p[1] 
err = p[2] 
bestfit = p[3]
phase = p[0] 
t = p[8]
ind = phase > 1.0
phase[ind] -= 1.0
resid =  np.zeros(len(phasebins)-1)

for j in range(1, len(phasebins)):
	#calculates uncertainty for in-eclipse points
	ind1 = (phase>=0.45)&(phase<=0.55)
	sig1 = np.sqrt(np.sum(err[ind1]**2)/sum(ind1)**2)

	#calculates uncertainty for phase bin
	ind2 = (phase >= phasebins[j-1]) & (phase < phasebins[j])
	sig2 = np.sqrt(np.sum(err[ind2]**2)/sum(ind2)**2)

	#x = np.mean(p[3][ind1])			#use for sine curve
	x = 1.						#use for physical model
	#print "make sure normalization is right"

	spectra[i, j-1, 0] = np.mean(data_corr[ind2])/x
	spectra[i, j-1, 1] = np.sqrt(sig1**2 + sig2**2)		#quadrature add uncertainties from baseline and bin 


#calculate beta factor to scale red noise
resid = data_corr - bestfit
period = 0.9255
binduration = (phasebins[2] - phasebins[1])*period
bins = np.arange(t.min(), t.max(), binduration)
binned_resid = np.zeros(len(bins)-1)
for ii in range(1,len(bins)-1):		
	ind = (t>bins[ii-1])&(t < bins[ii])
	binned_resid[ii] = np.mean(resid[ind]) 
	sigma = np.sqrt(np.sum(err[ind]**2)/sum(ind)**2)
print "sr, sw", np.sqrt(np.std(binned_resid[1::])**2 -  sigma**2), sigma
print "sr/sw", np.sqrt(np.std(binned_resid[1::])**2 -  sigma**2)/sigma
# s^2 = sw^2 + sr^2
sr = np.sqrt(np.std(binned_resid[1::])**2 -  sigma**2)

#adds red noise to bins
for j in range(1, len(phasebins)): spectra[i, j-1, 1] = np.sqrt(spectra[i, j-1,1]**2 + sr**2)

i += 1

#add Spitzer Ch 2
#f = "Ch2_best_fits/2017-10-11_20:24-zhang/bestfit.pic"
f = "Ch2_best_fits/2018-02-07_14:07-spherical/bestfit.pic"
waves.append(4.5)
dilution.append(0.1587)

p = pickle.load(open(f, 'rb'))
data_corr = p[1] 
err = p[2] 
bestfit = p[3]
phase = p[0] 
ind = phase > 1.0
phase[ind] -= 1.0

for j in range(1, len(phasebins)):
	#calculates uncertainty for in-eclipse points
	ind1 = (phase>=0.45)&(phase<=0.55)
	sig1 = np.sqrt(np.sum(err[ind1]**2)/sum(ind1)**2)

	#calculates uncertainty for phase bin
	ind2 = (phase >= phasebins[j-1]) & (phase < phasebins[j])
	sig2 = np.sqrt(np.sum(err[ind2]**2)/sum(ind2)**2)

#	x = np.mean(p[3][ind1])			#use for sine curve
	x = 1.						#use for physical model
	#print "make sure normalization is right"

	spectra[i, j-1, 0] = np.mean(data_corr[ind2])/x
	spectra[i, j-1, 1] = np.sqrt(sig1**2 + sig2**2)		#quadrature add uncertainties from baseline and bin 


i += 1
fig = plt.figure(figsize=(14,5))
ncol = 4
gs = gridspec.GridSpec(2, ncol, wspace = 0., hspace = 0.)

k = len(phasebins)
row = 0
counter = 1
for i in range(1, k):
	phase = (phasebins[i]+ phasebins[i-1])/2.

	if i == ncol+1: row += 1
	col = (i-1-2*row)%ncol

	ax = plt.subplot(gs[row, col])

	if (phase < 0.4)|(phase > 0.6):
		plt.errorbar(waves, (spectra[:,i-1,0]-1.)*1.0e3*(1+np.array(dilution)), yerr = 1.0e3*spectra[:,i-1,1], fmt='.k', zorder=100)	
		plt.gca().text(0.1, 0.8, '$\phi = $' + '{0:0.1f}'.format(phase), transform=ax.transAxes, fontsize=18)

		outname = "espec_phase_" + '{0:0.1f}'.format(phase) + ".txt"
		outfile = open(outname, "w")
		for j in range(len(dilution)):
			print>>outfile, waves[j], "{0:0.3e}".format((spectra[j, i-1,0]-1)*(1+dilution[j])), "{0:0.3e}".format(spectra[j,i-1,1]), phasebins[i-1], phasebins[i]
		outfile.close()


		rprs = 0.1146
		wave_hires, model_hires = best_fit_bb(waves, (spectra[:,i-1,0]-1.)*(1+np.array(dilution)), spectra[:,i-1,1], rprs)
		model_hires *= 1.e3
		plt.plot(wave_hires, model_hires, color='#6488ea', label='blackbody') 
                bbmin = np.min(model_hires)

                #plt.axhline(1, linestyle = 'dotted', color = '0.7')

                correction_factor = 1.096           #to account for the fact that Vivien used a different rp/rs
                if counter < 5: 
                    plt.plot(GCM[:,0], GCM[:,counter]*1e3*correction_factor, color = 'r')
                    print phase
                else: 
                    plt.plot(GCM[:,0], GCM[:,counter+1]*1e3*correction_factor, color = 'r', label = '$\\tau_\mathrm{drag4}$ GCM')
                    print phase
                counter += 1
	
                #LK
                if (col==0): plt.ylabel("$F_p$/$F_s$ ($\\times10^{-3}$)")
                if row ==1: plt.xlabel("Wavelength (microns)")

                if (col==3)&(row==1): plt.legend(loc=4, frameon=True) #legend

                plt.gca().set_yscale('log')
                plt.gca().set_xscale('log', basex=2)
                
                ax = plt.gca()
                ax.xaxis.set_major_locator(FixedLocator(np.array([1,2,4])))
                ax.xaxis.set_minor_locator(FixedLocator(np.array([1.1, 1.2, 1.3,  1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.2, 2.4, 2.6,  2.8, 3.0, 3.2,3.4, 3.6, 3.8, 4.4, 4.8])))
                ax.xaxis.set_minor_formatter(ticker.NullFormatter())
                ax.yaxis.set_minor_formatter(ticker.NullFormatter())
                ax.set_xticklabels(["1", "2", "4"])

                ax.yaxis.set_major_locator(FixedLocator(np.array([0.1, 0.2, 0.3, 0.4,  0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5])))
                ax.set_yticklabels(["0.1", "", "", "", "0.5", "", "", "", "", "1", "", "", "", "5"])

                plt.xlim(1.0, 5.0)
                
                """ymin = model_hires.min()*0.9 
                ymax = model_hires.max()*1.33
                if (col==3)&(row ==1): 
                        ymin *= 0.5
                        ymax = 2.5
                if (col==0)&(row ==0): ymin *= 0.7"""

                idx = row*4 + col
                ymin, ymax = ymins[idx], ymaxs[idx]
                plt.ylim(ymin, ymax)
                #ymin, ymax = plt.gca().get_ylim()
                yrange = np.log10(ymax) - np.log10(ymin)
                ystart = np.log10(ymin) + 0.1*yrange
                ystart = 10.**ystart
                print "row, col, yrange", row, col, np.log10(ymax)- np.log10(ymin)
                if (col==0)&(row==0): plt.gca().text(0.71, 0.13, '100\n ppm', transform=ax.transAxes, fontsize=14)
                plt.plot(np.array([4.5, 4.5]), np.array([ystart, ystart+ 0.1]), color = '0.5') 
                plt.plot(np.array([4.4, 4.6]), np.array([ystart, ystart]), color = '0.5') 
                plt.plot(np.array([4.4, 4.6]), np.array([ystart +0.1, ystart +0.1]), color = '0.5') 
                #plt.errorbar(4.2, ystart - 0.05, yerr = 0.05, color = '0.5', capsize = 10) 
                #plt.ylim(0.02, 5.1)

                if row==0: plt.gca().set_xticklabels([])
                
                
                #plt.savefig("spectrum"+ "{0:0.2f}".format(phase)+ ".png")
	#plt.show()

print "ERRORS & WARNINGS"
print "hard coding rprs for calculating best fit bb"
print "hard coding dilution in spitzer bandpass"
print "make sure normalization is right (sine curve vs. physical model)"
gs.tight_layout(fig, rect = [0., 0., 1., 1.])
plt.savefig("fig6.pdf")
