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

sns.set_context("notebook", font_scale = 1.3)
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
	#w, y, e = w[10], y[10], e[10]
#	w, y, e = w[11], y[11], e[11]

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
	outfile = open("temperatures.txt", "a")
	print>>outfile,  Tbest
	outfile.close()
	return waves_hires, blackbody(waves_hires*1.0e-6, Tbest)/star_bb_hires*rprs**2


path  = "WFC3_best_fits/spec_fits/"		#zhang model
files = glob.glob(os.path.join(path, "bestfit*.pic"))		

#phasebins = np.linspace(0.06, 0.44, 5)
#phasebins2 = np.linspace(0.56, 0.94, 5)

phasebins = np.array([0.06, 0.15, 0.25, 0.35, 0.44])
phasebins2 = np.array([0.56, 0.65, 0.75, 0.85, 0.94]) 

phasebins = np.append(phasebins, 0.5)
phasebins = np.append(phasebins, phasebins2)

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

		spectra[i, j-1, 0] = np.mean(data_corr[ind2])
		spectra[i, j-1, 1] = np.sqrt(sig1**2 + sig2**2)		#quadrature add uncertainties from baseline and bin 
		phases[j-1] = (phasebins[j] + phasebins[j-1])/2.

	i += 1

#add Spitzer Ch 1
f= "Ch1_best_fits/2017-10-11_20:25-zhang/bestfit.pic"
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
#print "sr, sw", np.sqrt(np.std(binned_resid[1::])**2 -  sigma**2), sigma
# s^2 = sw^2 + sr^2
sr = np.sqrt(np.std(binned_resid[1::])**2 -  sigma**2)

#adds red noise to bins
for j in range(1, len(phasebins)): spectra[i, j-1, 1] = np.sqrt(spectra[i, j-1,1]**2 + sr**2)

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
plt.figure(figsize=(14,5))
ncol = 4
gs = gridspec.GridSpec(2, ncol)

k = len(phasebins)
row = 0
for i in range(1, k):
	phase = (phasebins[i]+ phasebins[i-1])/2.

	if i == ncol+1: row += 1
	col = (i-1-2*row)%ncol

	ax = plt.subplot(gs[row, col])

	if (phase < 0.4)|(phase > 0.6):
		plt.errorbar(waves, (spectra[:,i-1,0]-1.)*1.0e3*(1+np.array(dilution)), yerr = 1.0e3*spectra[:,i-1,1], fmt='.k', zorder=100)	
		plt.gca().text(0.1, 0.8, '$\phi = $' + '{0:0.1f}'.format(phase), transform=ax.transAxes)

		outname = "espec_phase_" + '{0:0.1f}'.format(phase) + ".txt"
		outfile = open(outname, "w")
		for j in range(len(dilution)):
			print>>outfile, waves[j], "{0:0.3e}".format((spectra[j, i-1,0]-1)*(1+dilution[j])), "{0:0.3e}".format(spectra[j,i-1,1]), phasebins[i-1], phasebins[i]
		outfile.close()

		#w, f = rs.spectrum(phase, "all", "TiO-NoClouds.dat")
		#plt.plot(w, f*1.e3, color='#d3494e', linestyle='dashed', label= 'TiO-NoClouds')

		rprs = 0.115
		print "phase", phase
		wave_hires, model_hires = best_fit_bb(waves, (spectra[:,i-1,0]-1.)*(1+np.array(dilution)), spectra[:,i-1,1], rprs)
		model_hires *= 1.e3
		plt.plot(wave_hires, model_hires, color='#6488ea', label='blackbody') 
	

	if col==0: plt.ylabel("Planet-to-star flux (ppt)")
	if row ==1: plt.xlabel("Wavelength (microns)")

	if (col==3)&(row==1): plt.legend(loc=4)		#legend

	plt.gca().set_yscale('log')
	plt.gca().set_xscale('log', basex=2)
	
	ax = plt.gca()
	ax.xaxis.set_major_locator(FixedLocator(np.array([1,2,4])))
	ax.xaxis.set_minor_locator(FixedLocator(np.array([1.1, 1.2, 1.3,  1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.2, 2.4, 2.6,  2.8, 3.0, 3.2,3.4, 3.6, 3.8, 4.4, 4.8])))
	ax.xaxis.set_minor_formatter(ticker.NullFormatter())
	ax.yaxis.set_minor_formatter(ticker.NullFormatter())
	ax.set_xticklabels(["1", "2", "4"])


	ax.yaxis.set_major_locator(FixedLocator(np.array([1,2,3,4,5])))
	ax.set_yticklabels(["1", "2", "3", "4", "5"])


	plt.xlim(1.0, 5.0)
	ymin = model_hires.min()*0.9 
	ymax = model_hires.max()*1.3

	plt.ylim(ymin, ymax)
	if row==0: plt.gca().set_xticklabels([])
	
	
	#plt.savefig("spectrum"+ "{0:0.2f}".format(phase)+ ".png")
	#plt.show()

print "ERRORS & WARNINGS"
print "hard coding rprs for calculating best fit bb"
print "hard coding dilution in spitzer bandpass"
print "make sure normalization is right (sine curve vs. physical model)"
plt.savefig("emission_spectra.pdf")
