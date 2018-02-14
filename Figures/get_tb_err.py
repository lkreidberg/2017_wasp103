import pickle
import numpy as np
import matplotlib.pyplot as plt
import os, glob
from lc_fit import Model, LightCurveData
import matplotlib.gridspec as gridspec
from pylab import *
from matplotlib import rc, ticker
import scipy.stats as st
import emcee, corner

dayside = True

def quantile(x, q):
        return np.percentile(x, [100. * qi for qi in q])

#prior
def lnprior(theta):
        Ts, Tp = theta[0], theta[1]
	Ts_mu, Ts_err = 6110., 160.
	
	lnprior_prob = 0.
	if np.logical_or(Tp < 500., Tp > 4000.): lnprior_prob += - np.inf
        lnprior_prob -= 0.5*(((Ts - Ts_mu)/Ts_err)**2 + np.log(2.0*np.pi*(Ts_err)**2))
	return lnprior_prob
        
#likelihood function 
def lnlike(theta, x, y, err):
        Ts, Tp = theta[0], theta[1]
	rprs = 0.1146

	model = blackbody(x*1.e-6, Tp)/blackbody(x*1.e-6, Ts)*rprs**2 
        residuals = y - model

	#print 'Ts, Tp', Ts, Tp
	"""plt.errorbar(x,y,err, fmt='.k')
	plt.plot(x, model)
	plt.show()"""

        ln_likelihood = -0.5*(np.sum((residuals/err)**2 + np.log(2.0*np.pi*(err)**2)))

        return ln_likelihood

#posterior probability
def lnprob(theta, x, y, err):
        lp = lnprior(theta)
        if not np.isfinite(lp):
                return -np.inf
        return lp + lnlike(theta, x, y, err) 


def blackbody(l,T): 
	h = 6.62607e-34   #J/s
        c = 2.997925e8       #m/s
        k = 1.38065e-23     #J/K
        return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T))-1))                 #[W/sr/m^3]

path  = "WFC3_best_fits/spec_fits/"		#zhang model
files = glob.glob(os.path.join(path, "bestfit*.pic"))		

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
# s^2 = sw^2 + sr^2
sr = np.sqrt(np.std(binned_resid[1::])**2 -  sigma**2)

#adds red noise to bins
for j in range(1, len(phasebins)): spectra[i, j-1, 1] = np.sqrt(spectra[i, j-1,1]**2 + sr**2)

i += 1

#add Spitzer Ch 2
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

if dayside == True: 
	x, y, err,  best_fit_binned,  wl_hi,  y_hi_best, spec_arr, Tarr, P, samples = pickle.load(open("Mike_models/WASP-103b_grid_DAYSIDE_output.pic"))
	#x, y, err = x[0:10], y[0:10], err[0:10]
	x, y, err = x[10], y[10], err[10]
	#x, y, err = x[11], y[11], err[11]

	print x, y, err

	#run mcmc
	guess_Ts, guess_Tp = 6110., 1800.
	theta = [guess_Ts, guess_Tp]

	#initialize sampler
	ndim, nwalkers = len(theta), 50
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (x, y, err))
	pos = [theta + 1e-5*np.random.randn(ndim) for i in range(nwalkers)]

	#run mcmc
	sampler.run_mcmc(pos,2000)

	#plot MCMC output
	samples = sampler.chain[:, 1000:, :].reshape((-1, ndim))
	qs = quantile(samples[:,1], np.array([0.16, 0.5, 0.84]))
	print 0.5, qs[2] - qs[1] 
	#corner.corner(samples)
	#plt.show()


else:
	k = len(phasebins)
	for i in range(1, k):
		phase = (phasebins[i]+ phasebins[i-1])/2.

		if (phase < 0.4)|(phase > 0.6):
			x, y, err = np.array(waves), (spectra[:,i-1,0]-1.)*(1+np.array(dilution)), spectra[:,i-1,1]
			#x, y, err = x[0:10], y[0:10], err[0:10]
			#x, y, err = x[10], y[10], err[10]
			x, y, err = x[11], y[11], err[11]

			#run mcmc
			guess_Ts, guess_Tp = 6110., 1800.
			theta = [guess_Ts, guess_Tp]
		
			#initialize sampler
			ndim, nwalkers = len(theta), 50
			sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (x, y, err))
			pos = [theta + 1e-5*np.random.randn(ndim) for i in range(nwalkers)]

			#run mcmc
			sampler.run_mcmc(pos,2000)

			#plot MCMC output
			samples = sampler.chain[:, 1000:, :].reshape((-1, ndim))
			qs = quantile(samples[:,1], np.array([0.16, 0.5, 0.84]))
			print phase, qs[2] - qs[1] 
			#corner.corner(samples)
			#plt.show()

