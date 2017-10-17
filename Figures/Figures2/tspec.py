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

def blackbody(l,T):
        h = 6.626e-34   #J/s
        c = 3.0e8       #m/s
        kb = 1.38e-23     #J/K
        return 2*h*c**2/(l**5*(np.exp(h*c/(l*kb*T))-1))

def quantile(x, q):
        #return np.percentile(x, [100. * qi for qi in q])
        return np.percentile(x, 100. * q)

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

model = "PHYSICAL"
path  = "../" + model + "_MODEL_MCMC/"
files = glob.glob(os.path.join(path, "bestfit*.pic"))		

waves, dilution = [], []
rprs = np.zeros(len(files)+2)
rp_err = np.zeros(len(files)+2)
nightside = np.zeros(len(files)+2)

#WFC3
for i, f in enumerate(files):
	p = pickle.load(open(f, 'rb'))
	d, m, par = p[0], p[1], p[2]		#stores data,  model, and best fit parameters into d, m, & par

	mcmc_file = f[0:23] + "mcmc_out_" + f[31:35] + ".npy"
	mcmc = np.load(mcmc_file)

	phase = m.phase	
	waves.append(d.wavelength)
	dilution.append(d.dilution)
	
	nvisit = 4
	T_s, xi, T_n, delta_T, per, t0, eclipse =  par[d.par_order['T_s']*nvisit], par[d.par_order['xi']*nvisit], par[d.par_order['T_n']*nvisit], \
		par[d.par_order['delta_T']*nvisit], par[d.par_order['per']*nvisit], par[d.par_order['t0']*nvisit], False

	bestfit = np.array(spiderman_lc.lc(d.time, 0.115, T_s, d.l1, d.l2, xi, T_n, delta_T, False))
	
	ind1 = (phase>=0.95)|(phase<=0.05)
	nightside[i] = np.mean(bestfit[ind1]) - 1.

	rp = mcmc[:,0]

	rp_lo, rp_med, rp_hi = quantile(rp, np.array([0.16, 0.5, 0.84]))
	#rprs[i] = rp_med**2
	#rp_err[i] = 2. * rp_med * (((rp_med - rp_lo) + (rp_hi - rp_med))/2.) 
	rprs[i] = rp_med
	rp_err[i] = ((rp_med - rp_lo) + (rp_hi - rp_med))/2.

	#rprs[i] -= nightside[i]

	print d.wavelength, rprs[i], rp_err[i]


path  = "../TRANSIT_MCMC/"
files = glob.glob(os.path.join(path, "mcmc*.npy"))		

rprs2 = np.zeros(len(files)+2)
rp_err2 = np.zeros(len(files)+2)

#WFC3
for i, f in enumerate(files):
	mcmc = np.load(f)

	rp = mcmc[:,0]

	rp_lo, rp_med, rp_hi = quantile(rp, np.array([0.16, 0.5, 0.84]))
	rprs2[i] = rp_med**2
	rp_err2[i] = 2. * rp_med * (((rp_med - rp_lo) + (rp_hi - rp_med))/2.) 
	#rprs2[i] = rp_med
	#rp_err2[i] =   ((rp_med - rp_lo) + (rp_hi - rp_med))/2.

	#rprs[i] -= nightside[i]
	Tp = 800
	nightside = blackbody(waves[i]*1.e-6, Tp)/blackbody(waves[i]*1.e-6, 6110)*0.115

	print waves[i], rprs2[i], rp_err2[i], nightside

nightside = blackbody(4.5*1.e-6, Tp)/blackbody(4.5*1.e-6, 6110)*0.115
print 4.5, 0.12**2, 2*0.12*0.00128699949778, nightside


