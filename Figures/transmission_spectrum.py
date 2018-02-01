import numpy as np
import matplotlib.pyplot as plt
import os, glob
from lc_fit_nospiderman import Model, LightCurveData
import pickle
import pysynphot as psyn					#need to source activate astroconda
from matplotlib import rc, ticker
from pylab import *

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

path  = "WFC3_best_fits/spec_fits/"		#zhang model
files = glob.glob(os.path.join(path, "bestfit*.pic"))		
dilution = []
for f in files:
	p = pickle.load(open(f, 'rb'))
	d, m = p[0], p[1]			#stores data and model into d & m
	dilution.append(d.dilution)

dilution = np.array(dilution)

d = np.genfromtxt("w103b_transmission_spectrum.txt")

temps = [1400, 2000, 0]
labels = ["$T_\mathrm{night} = 1400$ K", "$T_\mathrm{night} = 2000$ K", "no correction"]
colors = ['blue', 'red', 'gray']

fig = plt.figure(figsize= (4,3))

for i, Tp in enumerate(temps):
	bb_planet, bb_star = psyn.BlackBody(Tp), psyn.BlackBody(6110.)

	nightside = bb_planet.flux/bb_star.flux*0.011 
	nightside_bin = np.interp(d[:,0]*1e4, bb_planet.wave, nightside) 
	
	if Tp == 0: nightside_bin = 0

	if i == 1: mean =  np.mean(d[:,1]*(1. + dilution)-nightside_bin)
	plt.errorbar(d[:,0], d[:,1]*(1. + dilution)-nightside_bin, yerr = d[:,2], marker = '.', linestyle = 'none', color = colors[i], label = labels[i], alpha = 0.7)

	# Spitzer Ch 2
	sp, sp_err = 0.11114407162, 0.00118749161363
	sp_err = sp*sp_err*2.
	#print "sp_err", sp_err*1e6
	sp = sp**2
	sp_dilution = 0.16
	nightside_bin = np.interp(4.5*1e4, bb_planet.wave, nightside) 
	if Tp == 0: nightside_bin = 0
	#plt.errorbar(4.5, sp*(1+sp_dilution), sp_err, fmt = '.r')
	plt.errorbar(4.5, sp*(1+sp_dilution) - nightside_bin, yerr = sp_err, xerr = 0.5, linestyle = 'none', color = colors[i], marker = '.', alpha = 0.7)

plt.gca().set_xscale('log')

ax = plt.gca()
ax.xaxis.set_major_locator(FixedLocator(np.array([1,2,4])))
ax.xaxis.set_minor_locator(FixedLocator(np.array([1.1, 1.2, 1.3,  1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.2, 2.4, 2.6,  2.8, 3.0, 3.2,3.4, 3.6, 3.8, 4.4, 4.8])))
ax.xaxis.set_minor_formatter(ticker.NullFormatter())
ax.yaxis.set_minor_formatter(ticker.NullFormatter())
ax.set_xticklabels(["1", "2", "4"])

plt.legend(loc = 'upper left')
plt.xlabel("Wavelength (microns)")
plt.ylabel("Transit Depth")
plt.axhline(mean, linestyle = 'dotted', zorder = -10)
ymin, ymax= 0.0128, 0.0147
plt.ylim(ymin, ymax)

#scale height
#kt/(mu*g) = 5.5e6 m (assuming mu = 2.3 amu, g = 15.85 m/s^2, T = 2410 K (phase 0.25 temperature)
# input to wolframalpha.com: 2*boltzmann constant*2410 Kelvin/(2.3 atomic mass units *15.85 m/s^2)*1.53 jupiter radii/(1.44 solar radii)^2
# feature amplitude = 1.2e-4

x2 = plt.gca().twinx()
scaleheight = 1.2e-4 
plt.plot(d[:,0], np.linspace(0, ymax - ymin, len(d[:,0]))/scaleheight, linewidth=0.)
plt.ylabel("scale heights")

plt.tight_layout()
plt.savefig("w103b_transmission_spectrum.pdf")
plt.show()
