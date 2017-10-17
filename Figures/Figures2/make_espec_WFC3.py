#import matplotlib
#matplotlib.use('ps')
import pickle
import numpy as np
import matplotlib.pyplot as plt
import os, glob
from lc_fit import Model, LightCurveData
import read_spectra_from_Vivien as rs
import matplotlib.gridspec as gridspec
from pylab import *
from matplotlib import rc

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

def blackbody(l,T):
        h = 6.62607e-34   #J/s
        c = 2.997925e8       #m/s
        k = 1.38065e-23     #J/K
        return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T))-1))

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
	waves_hires = np.linspace(1.1, 1.7, 100)
	star_bb_hires = np.interp(waves_hires, star[:,0], star[:,1])*1.e24/(waves_hires*np.pi*4.)
	print "best fit temp", Tbest
	return waves_hires, blackbody(waves_hires*1.0e-6, Tbest)/star_bb_hires*rprs**2


path  = "../PHYSICAL_MODEL_MCMC/"
#path  = "../PHYSICAL_MODEL_LSQ_6bins/"
#path  = "../SINE_MODEL_LSQ_6bins/"
files = glob.glob(os.path.join(path, "bestfit*.pic"))		

#phasebins = np.array([0.06, 0.15, 0.25, 0.35, 0.44])
#phasebins2 = np.array([0.56, 0.65, 0.75, 0.85, 0.94]) 

#this binning gives you quadrature (phase 0.25 and 0.75)
phasebins = np.array([0.06, 0.15, 0.2, 0.3, 0.44])
phasebins2 = np.array([0.56, 0.65, 0.7, 0.8, 0.94]) 

phasebins = np.append(phasebins, 0.5)
phasebins = np.append(phasebins, phasebins2)

spectra = np.zeros((len(files), len(phasebins)-1, 2))
phases = np.zeros(len(phasebins) -1 )
waves, dilution = [], []

for i, f in enumerate(files):
	p = pickle.load(open(f, 'rb'))
	d, m = p[0], p[1]			#stores data and model into d & m
	
	waves.append(d.wavelength)
	dilution.append(d.dilution)
	#print "dilution", d.dilution

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

		#print sum(ind2), spectra[i, j-1,0], spectra[i, j-1, 1]
		spectra[i, j-1, 0] = np.mean(data_corr[ind2])
		spectra[i, j-1, 1] = np.sqrt(sig1**2 + sig2**2)		#quadrature add uncertainties from baseline and bin 
		phases[j-1] = (phasebins[j] + phasebins[j-1])/2.

#print spectra
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
		plt.gca().text(0.1, 0.8, '$\phi = $' + '{0:0.2f}'.format(phase), transform=ax.transAxes)

		outname = "phase_" + '{0:0.2f}'.format(phase) + "_wfc3.txt"
		outfile = open(outname, "w")
		for j in range(len(dilution)):
			print>>outfile, waves[j], "{0:0.3e}".format((spectra[j, i-1,0]-1)*(1+dilution[j])), "{0:0.3e}".format(spectra[j,i-1,1]), phasebins[i-1], phasebins[i]

		w, f = rs.spectrum(phase, "all", "NoTiO-NoClouds.dat")
		plt.plot(w, f*1.e3, color='#d3494e', linestyle='dashed', label= 'NoTiO-NoClouds')

		rprs = 0.115
		print "phase", phase
		wave_hires, model_hires = best_fit_bb(waves, (spectra[:,i-1,0]-1.)*(1+np.array(dilution)), spectra[:,i-1,1], rprs)
		model_hires *= 1.e3
		plt.plot(wave_hires, model_hires, color='#6488ea', label='blackbody') 
	

	if col==0: plt.ylabel("Planet-to-star flux (ppt)")
	if row ==1: plt.xlabel("Wavelength (microns)")

	if (col==3)&(row==1): plt.legend(loc=4, fontsize=11)		#legend

	
	plt.xlim(1.1, 1.7)
	plt.ylim(0.8*np.min(model_hires), 1.2*np.max(model_hires))
	if row==0: plt.gca().set_xticklabels([])
	
	
	#plt.savefig("spectrum"+ "{0:0.2f}".format(phase)+ ".png")
	#plt.show()

print "ERRORS & WARNINGS"
print "hard coding rprs for calculating best fit bb"
print "make sure normalization is right (sine curve vs. physical model)"
plt.savefig("emission_spectra_wfc3.png")
plt.show()
