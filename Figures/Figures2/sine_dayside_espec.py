#import matplotlib
#matplotlib.use('ps')
import numpy as np
import pickle
import matplotlib.pyplot as plt
import os, glob
from lc_fit import Model, LightCurveData
import read_spectra_from_Vivien as rs
import matplotlib.gridspec as gridspec
from pylab import *
from matplotlib import rc
import spiderman_lc

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def blackbody(l,T): 
	h = 6.626e-34   #J/s
	c = 3.0e8       #m/s
	kb = 1.38e-23     #J/K
	return 2*h*c**2/(l**5*(np.exp(h*c/(l*kb*T))-1))


def best_fit_bb(w, y, e, rprs):
	Ts = np.linspace(300, 3600, 300)
	chibest = 10000.
	Tbest = 0.	
	Tstar = 6110.
	w = np.array(w)

	#get stellar spectrum
	star = np.genfromtxt("wasp103_sed_fluxes.out")
	star_bb = np.interp(w, star[:,0], star[:,1])*1.e24/(w*np.pi*4.)

	for T in Ts:
		model = blackbody(w*1.0e-6, T)/star_bb*rprs**2
		chi2 = np.sum((y - model)**2/e**2)
		if chi2 < chibest: 
			chibest, Tbest = chi2, T
	waves_hires = np.linspace(1.0, 5.0, 100)
	star_bb_hires = np.interp(waves_hires, star[:,0], star[:,1])*1.e24/(waves_hires*np.pi*4.)
	print "Best fit blackbody temp and chisq: ", Tbest, chibest/(len(e)-1)
	return waves_hires, blackbody(waves_hires*1.0e-6, Tbest)/star_bb_hires*rprs**2


#model = "PHYSICAL"
model = "SINE"
path  = "../" + model + "_MODEL_LSQ/"
files = glob.glob(os.path.join(path, "bestfit*.pic"))		

waves, dilution = [], []
fpfs = np.zeros(len(files)+2)
fp_err = np.zeros(len(files)+2)

for i, f in enumerate(files):
	p = pickle.load(open(f, 'rb'))
	d, m, par = p[0], p[1], p[2]		#stores data,  model, and best fit parameters into d, m, & par
	
	waves.append(d.wavelength)
	dilution.append(d.dilution)

	ind = d.err < 9.0e7			#indices of outliers
	
	err = d.err[ind]/d.flux[ind]		#normalized per point uncertainty 
	phase = m.phase[ind]
	data_corr = m.data_corr[ind]
	
	#calculates uncertainty for in-eclipse points
	ind1 = (phase>=0.46)&(phase<=0.55)
	sig1 = np.sqrt(np.sum(err[ind1]**2)/sum(ind1)**2)
	sig2 = np.sqrt(np.sum(err[~ind1]**2)/sum(~ind1)**2)

	fpfs[i] = np.mean(m.bestfit_no_eclipse[ind1])
	fp_err[i] = np.sqrt(sig1**2 + sig2**2)

i+= 1
waves.append(3.6)
dilution.append(0.1712)
#add Spitzer Ch 1
"""f = path + "ch1_bestfit.pic" 

p = pickle.load(open(f, 'rb'))
phase = p[0] 
data_corr = p[1] 
err = p[2] 
t = p[8]
bestpar = p[9]
ind = phase > 1.0
phase[ind] -= 1.0

T_s, xi, T_n, delta_T = bestpar[21], bestpar[24], bestpar[25], bestpar[26]
eclipse = False
bestfit = np.array(spiderman_lc.lc(t, 0.115, T_s, 3.e-6, 4.e-6, xi, T_n, delta_T, eclipse))

#calculates uncertainty for in-eclipse points
ind1 = (phase>=0.46)&(phase<=0.55)
sig1 = np.sqrt(np.sum(err[ind1]**2)/sum(ind1)**2)
sig2 = np.sqrt(np.sum(err[~ind1]**2)/sum(~ind1)**2)

if model == "SINE": x = np.mean(p[3][ind1])			#use for sine curve
else: x = 1.						#use for physical model

fpfs[i] = np.mean(bestfit[ind1]) 
fp_err[i] = np.sqrt(sig1**2 + sig2**2)

#calculate beta factor to scale red noise
resid = data_corr - p[3]
binduration = 0.1*per
bins = np.arange(t.min(), t.max(), binduration)
binned_resid = np.zeros(len(bins)-1)
for ii in range(1,len(bins)-1):		
	ind = (t>bins[ii-1])&(t < bins[ii])
	binned_resid[ii] = np.mean(resid[ind]) 
	sigma = np.sqrt(np.sum(err[ind]**2)/sum(ind)**2)
sr = np.sqrt(np.std(binned_resid[1::])**2 -  sigma**2)
#print "beta", sr/sigma

#adds red noise to bins
fp_err[i] = np.sqrt(fp_err[i]**2 + sr**2)"""

i += 1

#add Spitzer Ch 2
f = "/home/kreidberg/Projects/Data_reduction/Spitzer/WASP103_11099/Ch2/analysis/run/fgc/ap2750715/fit_phase_curve/2017-01-03_00:12-sine_model/bestfit.pic"  
waves.append(4.5)
dilution.append(0.1587)

p = pickle.load(open(f, 'rb'))
data_corr = p[1] 
err = p[2] 
bestfit = p[3]
phase = p[0] 
ind = phase > 1.0
phase[ind] -= 1.0
t = p[8]
bestpar = p[9]

#calculates uncertainty for in-eclipse points
ind1 = (phase>=0.46)&(phase<=0.55)
sig1 = np.sqrt(np.sum(err[ind1]**2)/sum(ind1)**2)
sig2 = np.sqrt(np.sum(err[~ind1]**2)/sum(~ind1)**2)

if model == "SINE": x = np.mean(p[3][ind1])			#use for sine curve

depth1 = bestpar[2]
depth2 = bestpar[8]

fpfs[i] = (depth1+depth2)/2. + 1.

"""plt.plot(phase, bestfit/x - 1., ',k')
plt.axhline(fpfs[i])
plt.axhline(fpfs[i]-1.e-4)
plt.axhline(fpfs[i]+ 1.e-4)
plt.show()"""

#fpfs[i] = np.mean(bestfit[ind1])/x

fp_err[i] = np.sqrt(sig1**2 + sig2**2)

######################################################################################
# PLOT
for i in range(len(waves)): print waves[i], (fpfs[i] - 1.)*(1.+np.array(dilution)[i]), fp_err[i]
plt.errorbar(waves, (fpfs-1.)*1e3*(1.+np.array(dilution)), yerr = fp_err*1.e3, fmt = '.k')

rprs = 0.115
#fp_err[5] = 1000.
#fp_err[8] = 1000.
fp_err[-2] = 1000.
wave_hires, model_hires = best_fit_bb(waves, (fpfs-1.)*(1.+np.array(dilution)), fp_err, rprs)
model_hires *= 1.e3
plt.plot(wave_hires, model_hires, color='#6488ea', label='blackbody') 
	
plt.ylabel("Planet-to-star flux (ppt)")
plt.xlabel("Wavelength (microns)")

plt.gca().set_yscale('log')
plt.gca().set_xscale('log', basex=2)
	
ax = plt.gca()
ax.xaxis.set_major_locator(FixedLocator(np.array([1,2,4])))
ax.xaxis.set_minor_locator(FixedLocator(np.array([1.1, 1.2, 1.3,  1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.2, 2.4, 2.6,  2.8, 3.0, 3.2,3.4, 3.6, 3.8, 4.4, 4.8])))
ax.set_xticklabels(["1", "2", "4"])


ax.yaxis.set_major_locator(FixedLocator(np.array([1,2,3,4,5])))
ax.set_yticklabels(["1", "2", "3", "4", "5"])


plt.xlim(1.0, 5.0)

ymin, ymax = 1, 6
plt.ylim(ymin, ymax)

#d = np.genfromtxt("fit_2016_12_29_17:01.txt")
#plt.errorbar(waves[0:11], d[:,1]*1e3*(1.+np.array(dilution))[0:11], fp_err[0:11]*1e3, fmt = '.r')

print "ERRORS & WARNINGS"
print "hard coding rprs for calculating best fit bb"
print "hard coding dilution in spitzer bandpass"
plt.savefig("dayside_spectrum.png")
plt.show()
