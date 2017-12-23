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
from astropy.convolution import Gaussian1DKernel, convolve
g = Gaussian1DKernel(stddev=20)

sns.set_context("paper", font_scale=1.2)
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

	w, y, e, = w[0:10], y[0:10], e[0:10]

	#get stellar spectrum
	star = np.genfromtxt("wasp103_sed_fluxes.out")
	star_bb = np.interp(w, star[:,0], star[:,1])*1.e24/(w*np.pi*4.)
        #Tstar = 6110.
        #star_bb = blackbody(w*1.0e-6, Tstar)
	outliers = 0.

	for T in Ts:
		model = blackbody(w*1.0e-6, T)/star_bb*rprs**2
		chi2 = np.sum((y - model)**2/e**2)
		if chi2 < chibest: 
			chibest, Tbest, outliers = chi2, T, (y-model)/e
	waves_hires = np.linspace(1.0, 5.0, 100)
	star_bb_hires = np.interp(waves_hires, star[:,0], star[:,1])*1.e24/(waves_hires*np.pi*4.)
	#print "Best fit blackbody temp, chisq, and outliers: ", Tbest, chibest/(len(e)-1), outliers
	print "Best fit blackbody temp, chi2: ", Tbest, chibest/(len(y) - 1.)
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
	bestfit = np.array(spiderman_lc.lc(d.time, 0.1127, T_s, d.l1, d.l2, xi, T_n, delta_T, d.dilution, eclipse))
	bestfit = bestfit[ind]

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
bestfit = np.array(spiderman_lc.lc(t, 0.1127, T_s, 3.15e-6, 3.95e-6, xi, T_n, delta_T, dilution[i], eclipse))


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
print "beta", sr/sigma

#adds red noise to bins
print  "error increases by", np.sqrt(fp_err[i]**2 + sr**2)/fp_err[i]
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
bestfit = np.array(spiderman_lc.lc(t, 0.1127, T_s, 4.e-6, 5.e-6, xi, T_n, delta_T, dilution[i], eclipse))

fpfs[i] = np.mean(bestfit[ind1])


fp_err[i] = np.sqrt(sig1**2 + sig2**2)

######################################################################################
# PLOT
plt.figure(figsize = (4,3))

plt.errorbar(waves, (fpfs-1.)*1e3, yerr = fp_err*1.e3, fmt = 'ow', zorder=1000, ecolor = 'k', markeredgecolor = 'k', markeredgewidth = 1.0)

outfile = open("espec_dayside.txt", "w")
for i in range(len(waves)): print>>outfile, "{0:0.3f}".format(waves[i]), "\t", "{0:0.3e}".format(fpfs[i] - 1.), "{0:0.3e}".format(fp_err[i])
outfile.close()

rprs = 0.1127
wave_hires, model_hires = best_fit_bb(waves, fpfs-1., fp_err, rprs)
model_hires *= 1.e3
plt.plot(wave_hires, model_hires, color='0.5', label='blackbody', linestyle = 'dashed', zorder = -20)
	

#fits from Mike
#wl,y_low_2sig, y_low_1sig, y_median, y_high_1sig, y_high_2sig=pickle.load(open("Retrieval/WASP103b_DAYSIDE_NOMINAL_spec.pic", "rb"))
data_wl, data, data_err,  best_fit_binned,  wl_hi,  y_hi_best, spec_arr, Tarr, P, samples = pickle.load(open("Mike_models/WASP-103b_grid_DAYSIDE_output.pic"))

#plt.fill_between(wl[::-1], convolve(y_low_2sig, g), convolve(y_high_2sig,g), color = 'orange', alpha = 0.5, zorder = -11)
#plt.fill_between(wl[::-1], convolve(y_low_1sig, g), convolve(y_high_1sig,g), color = 'orange', alpha = 0.5, zorder = -10)
#plt.plot(wl[::-1], convolve(y_median, g), zorder = -9, label = 'best fit')

plt.plot(wl_hi, y_hi_best, zorder = -9, label = 'best fit')

#scale = 1.1 
#w, f = rs.spectrum(0.5, "all", "TiO-NoClouds-Drag1.dat")
#plt.plot(w, f*1.e3*scale, color='#d3494e', linestyle='dashed', label= 'GCM with TiO')

#ch2_band = rs.spectrum(0.5, 'Spitzer2', "NoTiO-NoClouds.dat")
#plt.errorbar(4.5, ch2_band*1.e03, xerr = [0.5], color = "#d3494e", marker = 's')

plt.legend(loc = 'lower right', frameon=True, fontsize = 11)

plt.ylabel("Planet-to-star flux (ppt)")
plt.xlabel("Wavelength (microns)")

plt.xlim(1.1, 1.7)
ymin, ymax = 1.2, 2
plt.ylim(ymin, ymax)

#plots Cartier spectrum
"""s = np.genfromtxt("cartier_dayside_spec.txt")
plt.errorbar(s[:,0], s[:,1]/1.e3, s[:,3]/1.e3, fmt = 'xk')

#plots LK spectrum
s = np.genfromtxt("kreidberg13360_dayside_spec.txt")
plt.errorbar(s[:,0], s[:,1], s[:,2], fmt = 'xr')"""

a = plt.axes([.23, .62, .2, .28]) 
wl,y_low_2sig, y_low_1sig, y_median, y_high_1sig, y_high_2sig=pickle.load(open("Retrieval/WASP103b_DAYSIDE_NOMINAL_spec.pic", "rb"))
#wl,y_low_2sig, y_low_1sig, y_median, y_high_1sig, y_high_2sig=pickle.load(open("Mike_models/WASP-103b_grid_DAYSIDE_output.pic", "rb"))
plt.fill_between(wl[::-1], convolve(y_low_2sig, g), convolve(y_high_2sig,g), color = 'orange', alpha = 0.5, zorder = -11)
plt.fill_between(wl[::-1], convolve(y_low_1sig, g), convolve(y_high_1sig,g), color = 'orange', alpha = 0.5, zorder = -10)
plt.plot(wl[::-1], convolve(y_median, g), zorder = -9, label = 'best fit')

plt.plot(wave_hires, model_hires, color='0.5', label='blackbody', linestyle = 'dashed', zorder = -20)

plt.errorbar(waves, (fpfs-1.)*1e3, yerr = fp_err*1.e3, fmt = 'ow', zorder=1000, ecolor = 'k', markeredgecolor = 'k', markeredgewidth = 1.0)
#print (fpfs-1.)*1.e3
#print fp_err*1e3

plt.xlim(3, 5)
plt.ylim(3.5,7)


print "NOTE:"
print "you are hard coding rprs for calculating best fit bb"
plt.tight_layout()
plt.savefig("dayside_spectrum.pdf")
plt.show()
