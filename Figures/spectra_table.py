import pickle
import numpy as np
import matplotlib.pyplot as plt
import os, glob
from lc_fit import Model, LightCurveData


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

n = len(dilution)

print "\\begin{deluxetable*}{llllllllllll}"
#\tabletypesize{\footnotesize} 
print "\\tablecolumns{12}"
print "\\tablewidth{0pt}:"
print "\\tablecaption{Phase-Resolved Emission Spectra \label{table:spectra}}"
print "\\tablehead{"
print "\colhead{$\lambda$} & \colhead{Dilution} & \colhead{0.1} & \colhead{0.2} & \colhead{0.3} & \colhead{0.4} & \colhead{0.6} & \colhead{0.7} & \colhead{0.8} & \colhead{0.9} & \colhead{X} & \colhead{Y}}"
print "\startdata"

for j in range(len(dilution)):
	print waves[j], "&", "{0:0.2f}".format(dilution[j]), "&",
	for i in range(1, len(phasebins)):
		if i < len(phasebins)-1: print "$", "{0:0d}".format(int((spectra[j, i-1,0]-1)*(1+dilution[j])*1e6)), "\\pm", "{0:0d}".format(int(spectra[j,i-1,1]*1e6)), "$", "&", 
		else: print "$", "{0:0d}".format(int((spectra[j, i-1,0]-1)*(1+dilution[j])*1e6)), "\\pm", "{0:0d}".format(int(spectra[j,i-1,1]*1e6)), "$", "\\\\",
		
	print ""

print "\enddata"
print "\\vspace{-0.8cm}"
print "\\tablecomments{comments}"
print "\end{deluxetable*}"



print "ERRORS & WARNINGS"
print "hard coding rprs for calculating best fit bb"
print "hard coding dilution in spitzer bandpass"
