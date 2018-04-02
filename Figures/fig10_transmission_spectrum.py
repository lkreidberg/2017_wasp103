import numpy as np
import matplotlib.pyplot as plt
import os, glob
from lc_fit_nospiderman import Model, LightCurveData
import pickle
import pysynphot as psyn					#need to source activate astroconda
from matplotlib import rc, ticker
from pylab import *
import scipy.stats as st

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def get_significance(chisq, dof):
        alpha = (1. - st.chi2.cdf(chisq, dof))/2.
        z = st.norm.ppf(1.-alpha)
        return z

def weighted_mean(data, err):                #calculates the weighted mean for data points data with variances err
        weights = 1.0/err
        mu = np.sum(data*weights)/np.sum(weights)
        var = 1.0/np.sum(weights)
        return [mu, var]                #returns weighted mean and variance

path  = "WFC3_best_fits/spec_fits/"		#zhang model
files = glob.glob(os.path.join(path, "bestfit*.pic"))		
dilution = []
for f in files:
	p = pickle.load(open(f, 'rb'))
	d, m = p[0], p[1]			#stores data and model into d & m
	dilution.append(d.dilution)

dilution = np.array(dilution)

d = np.genfromtxt("w103b_transmission_spectrum.txt")
d_nodilution = np.genfromtxt("w103b_transmission_spectrum_040218_nodilution.txt")
d[:,1] = d_nodilution[:,1]      #uses least-squares best fit transit depths from a fit that did NOT including dilution


#temps = [1400, 1800, 0]
#labels = ["$T_\mathrm{night} = 1400$ K", "$T_\mathrm{night} = 1800$ K", "no correction"]
temps = [1700, 0]
labels = ["$T_\mathrm{night} = 1700$ K", "no correction"]
colors = ['blue', 'red', 'gray']

fig = plt.figure(figsize= (5,3))
rprs2 = 0.01313

wave, tspec, tspec_err = np.zeros(12), np.zeros(12), np.zeros(12)

for i, Tp in enumerate(temps):
	bb_planet, bb_star = psyn.BlackBody(Tp), psyn.BlackBody(6110.)

	nightside = bb_planet.flux/bb_star.flux*rprs2
	nightside_bin = np.interp(d[:,0]*1e4, bb_planet.wave, nightside) 
	
	if Tp == 0: nightside_bin = 0

	if i == 1: mean =  np.mean(d[:,1]*(1. + dilution)-nightside_bin)
        print "mean WFC3 depth = ", mean
	plt.errorbar(d[:,0], d[:,1]*(1. + dilution)-nightside_bin, yerr = d[:,2], marker = '.', linestyle = 'none', color = colors[i], label = labels[i], alpha = 0.7)
        tspec[0:10] =  d[:,1]*(1. + dilution)-nightside_bin 
        tspec_err[0:10] = d[:,2]
	wave[0:10] = d[:,0]

	# Spitzer Ch 1
	sp, sp_err = 0.10944, 0.0015
	sp_err = sp*sp_err*2.
	sp = sp**2
	sp_dilution = 0.17
	nightside_bin = np.interp(3.6*1e4, bb_planet.wave, nightside) 
	if Tp == 0: nightside_bin = 0
	plt.errorbar(3.6, sp*(1+sp_dilution) - nightside_bin, yerr = sp_err, xerr = 0.5, linestyle = 'none', color = colors[i], marker = '.', alpha = 0.7)
	tspec[10] =  sp*(1+sp_dilution) - nightside_bin 
        tspec_err[10] = sp_err
	wave[10] = 3.6

	# Spitzer Ch 2
	sp, sp_err = 0.11114407162, 0.00118749161363
	sp_err = sp*sp_err*2.
	#print "sp_err", sp_err*1e6
	sp = sp**2
	sp_dilution = 0.16
	nightside_bin = np.interp(4.5*1e4, bb_planet.wave, nightside) 
	if Tp == 0: nightside_bin = 0
	plt.errorbar(4.5, sp*(1+sp_dilution) - nightside_bin, yerr = sp_err, xerr = 0.5, linestyle = 'none', color = colors[i], marker = '.', alpha = 0.7)
	tspec[11] =  sp*(1+sp_dilution) - nightside_bin 
        tspec_err[11] = sp_err
	wave[11] = 4.5

        #gets reduced chi-sq for straight line fit
        mu, sig = weighted_mean(tspec, tspec_err**2)
        chi2, dof  = np.sum((tspec-mu)**2/tspec_err**2), (len(tspec) - 1)
        print "temperature, flat line rejection confidence = ", Tp, get_significance(chi2, dof)
        print "rprs mean, err", mu, np.sqrt(sig)*1e2

	for i in range(len(tspec)):
		print wave[i], tspec[i], tspec_err[i]

#plot GCM
d = np.genfromtxt("W103b-Drag4-transmission-spectrum-CaTiO3-1microns.dat", skip_header = 1, delimiter = ',')
rprs = (d[:,1] + 1.135e8)/9.99e8
rprs2 = rprs**2
plt.plot(d[:,0], rprs2, color = 'k', zorder = -10, linewidth = 1.2, label = "$\\tau_\mathrm{drag4}$ GCM",alpha = 0.7) 

ind = (4. < d[:,0])&(d[:,0] < 5.)
plt.plot(4.5, np.mean(rprs2[ind]), marker = 's', color = 'k', alpha = 0.7)
ind = (3.1 < d[:,0])&(d[:,0] < 4.1)
plt.plot(3.6, np.mean(rprs2[ind]), marker = 's', color = 'k', alpha = 0.7)

"""cloudrprs = (d[:,2] + 1.125e8)/9.99e8
cloudrprs2 = cloudrprs**2
print np.mean(cloudrprs2)
plt.plot(d[:,0], cloudrprs2, color = '0.3', linestyle = 'dotted', zorder = -10)"""

plt.gca().set_xscale('log')

ax = plt.gca()
ax.xaxis.set_major_locator(FixedLocator(np.array([1,2,4])))
ax.xaxis.set_minor_locator(FixedLocator(np.array([1.1, 1.2, 1.3,  1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.2, 2.4, 2.6,  2.8, 3.0, 3.2,3.4, 3.6, 3.8, 4.4, 4.8])))
ax.xaxis.set_minor_formatter(ticker.NullFormatter())
ax.yaxis.set_minor_formatter(ticker.NullFormatter())
ax.set_xticklabels(["1", "2", "4"])
ax.yaxis.set_major_locator(FixedLocator(np.array([0.0130, 0.0135, 0.0140, 0.0145])))


plt.legend(loc = 'upper left')
plt.xlabel("Wavelength (microns)")
plt.ylabel("Transit Depth")
#plt.axhline(mean, linestyle = 'dotted', zorder = -10)
ymin, ymax= 0.0128, 0.0147
plt.ylim(ymin, ymax)

#scale height
#kt/(mu*g) = 5.5e6 m (assuming mu = 2.3 amu, g = 15.85 m/s^2, T = 2410 K (phase 0.25 temperature)
# input to wolframalpha.com: 2*boltzmann constant*2410 Kelvin/(2.3 atomic mass units *15.85 m/s^2)*1.53 jupiter radii/(1.44 solar radii)^2
# feature amplitude = 1.2e-4

x2 = plt.gca().twinx()
scaleheight = 1.2e-4 
plt.plot(d[:,0], np.linspace(0, ymax - ymin, len(d[:,0]))/scaleheight - (mean-ymin)/scaleheight, linewidth=0.)
plt.ylabel("Scale Height")

plt.xlim(1,5.2)
plt.tight_layout()
plt.savefig("fig10.pdf")
