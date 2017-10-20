import numpy as np
import pickle
import matplotlib.pyplot as plt
import os, glob
from lc_fit import Model, LightCurveData
#import read_spectra_from_Vivien as rs
import matplotlib.gridspec as gridspec
from pylab import *
from matplotlib import rc
import spiderman_lc
from matplotlib import ticker
import seaborn as sns

sns.set_context("talk", font_scale=1.0)
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})


rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def blackbody(l,T): 
	h = 6.62607e-34   #J/s
        c = 2.997925e8       #m/s
        k = 1.38065e-23     #J/K
        return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T))-1))                 #[W/sr/m^3]

def best_fit_bb(w, y, e, rprs):
	Ts = np.linspace(300, 3600, 300)
	chibest = 10000.
	Tbest = 0.	
	Tstar = 6110.
	w = np.array(w)

	#w, y, e = w[0:10], y[0:10], e[0:10]
	#w, y, e = w[11], y[11], e[11]
	#print "just fitting WFC3 data"

	#get stellar spectrum
	star = np.genfromtxt("wasp103_sed_fluxes.out")
	star_bb = np.interp(w, star[:,0], star[:,1])*1.e24/(w*np.pi*4.)
	outliers = 0.

	for T in Ts:
		model = blackbody(w*1.0e-6, T)/star_bb*rprs**2
		chi2 = np.sum((y - model)**2/e**2)
		if chi2 < chibest: 
			chibest, Tbest, outliers = chi2, T, (y-model)/e
	waves_hires = np.linspace(1.0, 5.0, 100)
	star_bb_hires = np.interp(waves_hires, star[:,0], star[:,1])*1.e24/(waves_hires*np.pi*4.)
	#print "Best fit blackbody temp, chisq, and outliers: ", Tbest, chibest/(len(e)-1), outliers
	print "Best fit blackbody temp: ", Tbest
	return waves_hires, blackbody(waves_hires*1.0e-6, Tbest)/star_bb_hires*rprs**2


model = "PHYSICAL"
path  = "WFC3_best_fits/spec_fits/"
files = glob.glob(os.path.join(path, "bestfit*.pic"))		

waves, dilution = [], []
nspitzchan = 2
fpfs = np.zeros(len(files)+nspitzchan)
fp_err = np.zeros(len(files)+nspitzchan)

for i, f in enumerate(files):
	p = pickle.load(open(f, 'rb'))
	d, m, par = p[0], p[1], p[2]		#stores data,  model, and best fit parameters into d, m, & par
	
	waves.append(d.wavelength)
	dilution.append(d.dilution)


	ind = d.err < 9.0e7			#indices of outliers
	
	err = d.err[ind]/d.flux[ind]		#normalized per point uncertainty 
	phase = m.phase[ind]
	data_corr = m.data_corr[ind]
	nvisit = 4
	T_s, xi, T_n, delta_T, per, t0, eclipse =  par[d.par_order['T_s']*nvisit], par[d.par_order['xi']*nvisit], par[d.par_order['T_n']*nvisit], \
		par[d.par_order['delta_T']*nvisit], par[d.par_order['per']*nvisit], par[d.par_order['t0']*nvisit], False
	bestfit = np.array(spiderman_lc.lc(d.time, 0.115, T_s, d.l1, d.l2, xi, T_n, delta_T, d.dilution, eclipse))
	bestfit = bestfit[ind]

	#print d.wavelength, T_n, delta_T

	#t_eclipse = np.linspace(0.45, 0.55, 100)*per + t0
	#bestfit_eclipse = np.array(spiderman_lc.lc(t_eclipse, 0.115, T_s, d.l1, d.l2, xi, T_n, delta_T, eclipse))
	#LK note: this only makes a difference of 3 ppm
	
	#calculates uncertainty for in-eclipse points
	ind1 = (phase>=0.46)&(phase<=0.55)
	sig1 = np.sqrt(np.sum(err[ind1]**2)/sum(ind1)**2)
	sig2 = np.sqrt(np.sum(err[~ind1]**2)/sum(~ind1)**2)

	fpfs[i] = np.mean(bestfit[ind1])
	fp_err[i] = np.sqrt(sig1**2 + sig2**2)

i+= 1
#add Spitzer Ch 1
f= "Ch1_best_fits/2017-10-11_20:25-zhang/bestfit.pic"
waves.append(3.6)
dilution.append(0.1712)

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
bestfit = np.array(spiderman_lc.lc(t, 0.115, T_s, 3.15e-6, 3.95e-6, xi, T_n, delta_T, dilution[i], eclipse))


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
fp_err[i] = np.sqrt(fp_err[i]**2 + sr**2)

i += 1

#add Spitzer Ch 2
f = "Ch2_best_fits/2017-10-11_20:24-zhang/bestfit.pic"
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


T_s, xi, T_n, delta_T = bestpar[21], bestpar[24], bestpar[25], bestpar[26]			#for lin rmap
#T_s, xi, T_n, delta_T = bestpar[22], bestpar[25], bestpar[26], bestpar[27]			#for quad ramp
eclipse = False
bestfit = np.array(spiderman_lc.lc(t, 0.115, T_s, 4.e-6, 5.e-6, xi, T_n, delta_T, dilution[i], eclipse))

fpfs[i] = np.mean(bestfit[ind1])


fp_err[i] = np.sqrt(sig1**2 + sig2**2)

######################################################################################
# PLOT
plt.errorbar(waves, (fpfs-1.)*1e3, yerr = fp_err*1.e3, fmt = '.k', zorder=1000)

outfile = open("espec_dayside.txt", "w")
for i in range(len(waves)): print>>outfile, "{0:0.3f}".format(waves[i]), "\t", "{0:0.3e}".format(fpfs[i] - 1.), "{0:0.3e}".format(fp_err[i])
outfile.close()

rprs = 0.115
#fp_err[5] = 1000.
#fp_err[8] = 1000.
wave_hires, model_hires = best_fit_bb(waves, fpfs-1., fp_err, rprs)
model_hires *= 1.e3
plt.plot(wave_hires, model_hires, color='#6488ea', label='blackbody') 
	

scale = 1.1 
#w, f = rs.spectrum(0.5, "all", "TiO-NoClouds-Drag1.dat")
#plt.plot(w, f*1.e3*scale, color='#d3494e', linestyle='dashed', label= 'GCM with TiO')

#ch2_band = rs.spectrum(0.5, 'Spitzer2', "NoTiO-NoClouds.dat")
#plt.errorbar(4.5, ch2_band*1.e03, xerr = [0.5], color = "#d3494e", marker = 's')

plt.legend(loc =2, frameon=True)

plt.ylabel("Planet-to-star flux (ppt)")
plt.xlabel("Wavelength (microns)")

plt.gca().set_yscale('log')
plt.gca().set_xscale('log', basex=2)
	
ax = plt.gca()
ax.xaxis.set_major_locator(FixedLocator(np.array([1,2,4])))
ax.xaxis.set_minor_locator(FixedLocator(np.array([1.1, 1.2, 1.3,  1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.2, 2.4, 2.6,  2.8, 3.0, 3.2,3.4, 3.6, 3.8, 4.4, 4.8])))
ax.set_xticklabels([1, 2, 4])
ax.xaxis.set_minor_formatter(ticker.NullFormatter())


ax.yaxis.set_major_locator(FixedLocator(np.array([1,2,3,4,5])))
ax.set_yticklabels(["1", "2", "3", "4", "5"])
ax.yaxis.set_minor_formatter(ticker.NullFormatter())


plt.xlim(1.0, 5.0)

ymin, ymax = 1, 6
plt.ylim(ymin, ymax)


print "ERRORS & WARNINGS"
print "hard coding rprs for calculating best fit bb - is this right?"
print "hard coding dilution in spitzer bandpass - are we also supposed to do this for WFC3?"
plt.savefig("dayside_spectrum.pdf")
#plt.show()
