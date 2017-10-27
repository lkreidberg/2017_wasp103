import sys
#sys.path.insert(1, '/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/')
import mpfit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import glob
import os
from astropy.io import ascii
import getopt
import seaborn as sns
import matplotlib
import matplotlib.cm as cm
from datetime import datetime
import time as pythontime
import emcee
#import corner
import batman
import pickle
import spiderman

sns.set_context("talk")
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})
matplotlib.rcParams.update({'lines.markeredgewidth':0.3})
matplotlib.rcParams.update({'axes.formatter.useoffset':False})

def quantile(x, q):
        return np.percentile(x, [100. * qi for qi in q])

def parse_skip_orbs(x):
	n = len(x)
	i = 0
	temp = 0
	skip_orbs = []
	while(i < n):
		if x[i] == "_": 
			skip_orbs.append(int(x[temp:i]))
			temp = i+1
		i += 1
	return np.array(skip_orbs)

class LightCurveData:
	"""
	doc
	"""
	def __init__(self, data_file, obs_par, fit_par, flags):
		par_order = {line['parameter']: i for i, line in enumerate(fit_par)}

		d = np.genfromtxt(data_file)

		ind = np.arange(len(d))
		d = d[(ind!=410)&(ind!=375)&(ind!=376)&(ind!=377)&(ind!=520)&(ind!=378)&(ind!=379)]
		ind = d[:,1]<0.5*np.median(d[:,1])
		d = d[~ind]			
		#print "LK quick fix: remove data point 410 (weird outlier) and points with flux < 0.5*median flux to get rid of bad points in 14th orbit. FIXME for other data sets"

		d = d[np.argsort(d[:,5])]				#sorts data by time

		use_first_exp = False
		if use_first_exp == False:
			ind = np.diff(d[:,5]) < 30./60./24.		#removes first exposure from each orbit
			d = d[1:][ind]

		orb_num = d[:,7]					#removes first orbit from each visit 

		if obs_par['lc_type'] == "transit": skip_orbs = parse_skip_orbs(obs_par['skip_orbs_transit'])
		elif obs_par['lc_type'] == "eclipse": skip_orbs = parse_skip_orbs(obs_par['skip_orbs_eclipse'])
		elif obs_par['lc_type'] == "phase_curve": skip_orbs = parse_skip_orbs(obs_par['skip_orbs_phase_curve'])
		elif obs_par['lc_type'] == "physical_model": skip_orbs = parse_skip_orbs(obs_par['skip_orbs_physical_model'])
		else: raise Exception("Unsupported lc_type")

		ind = orb_num>-1 
		for i in skip_orbs: ind[orb_num == i] = False 
		d = d[ind]
		orb_num = d[:,7]

		for i in range(len(orb_num)-1):				#compresses orb_num
			if orb_num[i+1] - orb_num[i] > 1: 
				ind = orb_num==orb_num[i+1]
				orb_num[ind] = orb_num[i] + 1
		orb_num -= orb_num[0]
		norbit = int(orb_num[-1]+1)
		
		n = len(d)
		vis_num = np.zeros(n)
		t_vis = np.zeros(n) 
		t_orb = np.zeros(n) 
		t_delay = np.zeros(n) 
		
		visit = 0
		for i in range(1,n):
			vis_num[i - 1] = visit
			if (d[i,5] - d[i-1,5]) > 9./24.: visit += 1	#stores data as separate visit if gap is longer than 9 hrs; FIXME: what if visits are closer together than 9 hours?
		vis_num[-1] = visit

		nvisit = int(obs_par['nvisit'])
		for i in range(nvisit - 2):
			ind = vis_num == i
			t_vis[ind] = d[ind,5] - d[ind,5][0]
		
		orbs_per_visit = norbit/nvisit

		for i in range(norbit):
			ind = orb_num == i
			t_orb[ind] = d[ind,5] - d[ind,5][0]
			if i%orbs_per_visit == False: t_delay[ind] = 1.

		err = np.sqrt(d[:,2])
		flux = d[:,1]
		time  = d[:,5]
		scan_direction = d[:,8]

		wavelength = d[0,3]

		#fixes limb darkening if "fix_ld" parameter is set to True in obs_par.txt
		if obs_par['fix_ld'].lower() == "true":
			ld = np.genfromtxt(obs_par['ld_file'])
			
			i = 0
			while(wavelength > ld[i,1]): i += 1
			
			u1 = ld[i, 3]
			u2 = ld[i, 4] 
			fit_par['value'][np.where(fit_par['parameter']=='u1')] = u1
			fit_par['fixed'][np.where(fit_par['parameter']=='u1')] = "true"
			fit_par['value'][np.where(fit_par['parameter']=='u2')] = u2
			fit_par['fixed'][np.where(fit_par['parameter']=='u2')] = "true"

		nfree_param = 0
		for i in range(len(fit_par)):
			if fit_par['fixed'][i].lower() == "false":
				if fit_par['tied'][i].lower() == "true": nfree_param += 1
				else: nfree_param += nvisit
		nfree_param -= 2	#subtracts 2 because we're fixing v2 at 0 for zhao data

		if flags['fit-white']: zhao_data = pickle.load(open("zhao_data_pickles_white/zhao_data_white.p", "rb"))
		else: zhao_data = pickle.load(open("zhao_data_pickles_6bins/zhao_data_"+"{0:0.3f}".format(wavelength)+".p", "rb"))
		print "choosing zhao pickles by hand; make sure you specified the right file"
		print "specifying prior on T_n and delta_T_cloud by hand - indices could be wrong later!"

		self.time = np.append(time, zhao_data.time)
		self.flux = np.append(flux, zhao_data.flux)
		self.err = np.append(err, zhao_data.err)
		self.wavelength = wavelength
		self.exp_time = float(obs_par['exp_time'])
		self.nvisit = nvisit 
		self.vis_num = np.append(vis_num, zhao_data.vis_num + 2)
		self.orb_num = np.append(orb_num, zhao_data.orb_num)
		self.scan_direction = np.append(scan_direction, zhao_data.scan_direction)
		self.t_vis = np.append(t_vis, zhao_data.t_vis)
		self.t_orb = np.append(t_orb, zhao_data.t_orb)
		self.t_delay = np.append(t_delay, zhao_data.t_delay)
		self.par_order = par_order
		self.nfree_param = nfree_param
		self.dof = n + zhao_data.npoints - nfree_param 
		self.npoints = n + zhao_data.npoints
		self.lc_type = obs_par['lc_type']
		if flags['fit-white']:
			self.l1 = (1.1)*1.0e-6		#FIXME
			self.l2 = (1.7)*1.0e-6		#FIXME
		else:
			self.l1 = (self.wavelength - 0.0225)*1.0e-6	#specific to 12 bins!
			self.l2 = (self.wavelength + 0.0225)*1.0e-6
		self.stellar_grid = spiderman.stellar_grid.gen_grid(self.l1, self.l2)
		print "specifying l1 and l2 by hand - FIXME"
		print "setting rp by hand in spiderman functions"
		self.all_sys = None
		self.prior = format_prior_for_mcmc(obs_par, fit_par)

		dil = np.genfromtxt("flux_contam_full.out", skip_header=2)
		ind = (dil[:,0] > self.l1*1.e6)&(dil[:,0]<self.l2*1.e6)

		self.dilution = np.mean(dil[ind,1])
		print "dilution", self.dilution


class FormatParams: 
	"""
	doc
	"""
	def __init__(self, params, data):
		self.per = params[data.par_order['per']*data.nvisit:(1 + data.par_order['per'])*data.nvisit]
		self.t0 = params[data.par_order['t0']*data.nvisit:(1 + data.par_order['t0'])*data.nvisit]
		self.t_secondary = params[data.par_order['t_secondary']*data.nvisit:(1 + data.par_order['t_secondary'])*data.nvisit]
		self.w = params[data.par_order['w']*data.nvisit:(1 + data.par_order['w'])*data.nvisit]
		self.a = params[data.par_order['a']*data.nvisit:(1 + data.par_order['a'])*data.nvisit]
		self.inc = params[data.par_order['inc']*data.nvisit:(1 + data.par_order['inc'])*data.nvisit]
		self.rp = params[data.par_order['rp']*data.nvisit:(1 + data.par_order['rp'])*data.nvisit]
		self.fp = params[data.par_order['fp']*data.nvisit:(1 + data.par_order['fp'])*data.nvisit]
		self.u1 = params[data.par_order['u1']*data.nvisit:(1 + data.par_order['u1'])*data.nvisit]
		self.u2 = params[data.par_order['u2']*data.nvisit:(1 + data.par_order['u2'])*data.nvisit]
		self.ecc = params[data.par_order['ecc']*data.nvisit:(1 + data.par_order['ecc'])*data.nvisit]
		self.c = params[data.par_order['c']*data.nvisit:(1 + data.par_order['c'])*data.nvisit]
		self.v = params[data.par_order['v']*data.nvisit:(1 + data.par_order['v'])*data.nvisit]
		self.v2 = params[data.par_order['v2']*data.nvisit:(1 + data.par_order['v2'])*data.nvisit]
		self.r1 = params[data.par_order['r1']*data.nvisit:(1 + data.par_order['r1'])*data.nvisit]
		self.r2 = params[data.par_order['r2']*data.nvisit:(1 + data.par_order['r2'])*data.nvisit]
		self.r3 = params[data.par_order['r3']*data.nvisit:(1 + data.par_order['r3'])*data.nvisit]
		self.scale = params[data.par_order['scale']*data.nvisit:(1 + data.par_order['scale'])*data.nvisit]
		self.amp1 = params[data.par_order['amp1']*data.nvisit:(1 + data.par_order['amp1'])*data.nvisit]
		self.theta1 = params[data.par_order['theta1']*data.nvisit:(1 + data.par_order['theta1'])*data.nvisit]
		self.amp2 = params[data.par_order['amp2']*data.nvisit:(1 + data.par_order['amp2'])*data.nvisit]
		self.theta2 = params[data.par_order['theta2']*data.nvisit:(1 + data.par_order['theta2'])*data.nvisit]
		self.a_abs = params[data.par_order['a_abs']*data.nvisit:(1 + data.par_order['a_abs'])*data.nvisit]
 		self.xi = params[data.par_order['xi']*data.nvisit:(1 + data.par_order['xi'])*data.nvisit]
 		self.T_n = params[data.par_order['T_n']*data.nvisit:(1 + data.par_order['T_n'])*data.nvisit]
 		self.delta_T = params[data.par_order['delta_T']*data.nvisit:(1 + data.par_order['delta_T'])*data.nvisit]
 		self.T_s = params[data.par_order['T_s']*data.nvisit:(1 + data.par_order['T_s'])*data.nvisit]
		self.p_u1 = params[data.par_order['p_u1']*data.nvisit:(1 + data.par_order['p_u1'])*data.nvisit]
		self.p_u2 = params[data.par_order['p_u2']*data.nvisit:(1 + data.par_order['p_u2'])*data.nvisit]
		self.delta_T_cloud = params[data.par_order['delta_T_cloud']*data.nvisit:(1 + data.par_order['delta_T_cloud'])*data.nvisit]
		self.la0 = params[data.par_order['la0']*data.nvisit:(1 + data.par_order['la0'])*data.nvisit]
		self.lo0 = params[data.par_order['lo0']*data.nvisit:(1 + data.par_order['lo0'])*data.nvisit]
		self.size = params[data.par_order['size']*data.nvisit:(1 + data.par_order['size'])*data.nvisit]
		self.spot_T = params[data.par_order['spot_T']*data.nvisit:(1 + data.par_order['spot_T'])*data.nvisit]
		self.p_T = params[data.par_order['p_T']*data.nvisit:(1 + data.par_order['p_T'])*data.nvisit]
		self.sph0 = params[data.par_order['sph0']*data.nvisit:(1 + data.par_order['sph0'])*data.nvisit]
		self.sph1 = params[data.par_order['sph1']*data.nvisit:(1 + data.par_order['sph1'])*data.nvisit]
		self.sph2 = params[data.par_order['sph2']*data.nvisit:(1 + data.par_order['sph2'])*data.nvisit]
		self.sph3 = params[data.par_order['sph3']*data.nvisit:(1 + data.par_order['sph3'])*data.nvisit]

def PrintParams(m, data): 
	print "per\t", m.params[data.par_order['per']*data.nvisit:(1 + data.par_order['per'])*data.nvisit]
	print "t0\t", m.params[data.par_order['t0']*data.nvisit:(1 + data.par_order['t0'])*data.nvisit]
	print "t_s\t", m.params[data.par_order['t_secondary']*data.nvisit:(1 + data.par_order['t_secondary'])*data.nvisit]
	print "w\t", m.params[data.par_order['w']*data.nvisit:(1 + data.par_order['w'])*data.nvisit]
	print "a\t", m.params[data.par_order['a']*data.nvisit:(1 + data.par_order['a'])*data.nvisit]
	print "inc\t", m.params[data.par_order['inc']*data.nvisit:(1 + data.par_order['inc'])*data.nvisit]
	print "rp\t", m.params[data.par_order['rp']*data.nvisit:(1 + data.par_order['rp'])*data.nvisit]
	print "fp\t", m.params[data.par_order['fp']*data.nvisit:(1 + data.par_order['fp'])*data.nvisit]
	print "u1\t", m.params[data.par_order['u1']*data.nvisit:(1 + data.par_order['u1'])*data.nvisit]
	print "u2\t", m.params[data.par_order['u2']*data.nvisit:(1 + data.par_order['u2'])*data.nvisit]
	print "ecc\t", m.params[data.par_order['ecc']*data.nvisit:(1 + data.par_order['ecc'])*data.nvisit]
	print "c\t", m.params[data.par_order['c']*data.nvisit:(1 + data.par_order['c'])*data.nvisit]
	print "v\t", m.params[data.par_order['v']*data.nvisit:(1 + data.par_order['v'])*data.nvisit]
	print "v2\t",  m.params[data.par_order['v2']*data.nvisit:(1 + data.par_order['v2'])*data.nvisit]
	print "r1\t", m.params[data.par_order['r1']*data.nvisit:(1 + data.par_order['r1'])*data.nvisit]
	print "r2\t", m.params[data.par_order['r2']*data.nvisit:(1 + data.par_order['r2'])*data.nvisit]
	print "r3\t", m.params[data.par_order['r3']*data.nvisit:(1 + data.par_order['r3'])*data.nvisit]
	print "scale\t", m.params[data.par_order['scale']*data.nvisit:(1 + data.par_order['scale'])*data.nvisit]
	print "amp1\t", m.params[data.par_order['amp1']*data.nvisit:(1 + data.par_order['amp1'])*data.nvisit]
	print "theta1\t", m.params[data.par_order['theta1']*data.nvisit:(1 + data.par_order['theta1'])*data.nvisit]
	print "amp2\t", m.params[data.par_order['amp2']*data.nvisit:(1 + data.par_order['amp2'])*data.nvisit]
	print "theta2\t", m.params[data.par_order['theta2']*data.nvisit:(1 + data.par_order['theta2'])*data.nvisit]
	print "a_abs\t", m.params[data.par_order['a_abs']*data.nvisit:(1 + data.par_order['a_abs'])*data.nvisit]
 	print "xi\t", m.params[data.par_order['xi']*data.nvisit:(1 + data.par_order['xi'])*data.nvisit]
 	print "T_n\t", m.params[data.par_order['T_n']*data.nvisit:(1 + data.par_order['T_n'])*data.nvisit]
 	print "delta_T\t", m.params[data.par_order['delta_T']*data.nvisit:(1 + data.par_order['delta_T'])*data.nvisit]
 	print "T_s\t", m.params[data.par_order['T_s']*data.nvisit:(1 + data.par_order['T_s'])*data.nvisit]
 	print "p_u1\t", m.params[data.par_order['p_u1']*data.nvisit:(1 + data.par_order['p_u1'])*data.nvisit]
 	print "p_u2\t", m.params[data.par_order['p_u2']*data.nvisit:(1 + data.par_order['p_u2'])*data.nvisit]
 	print "delta_T_cloud\t", m.params[data.par_order['delta_T_cloud']*data.nvisit:(1 + data.par_order['delta_T_cloud'])*data.nvisit]
 	print "la0\t", m.params[data.par_order['la0']*data.nvisit:(1 + data.par_order['la0'])*data.nvisit]
 	print "lo0\t", m.params[data.par_order['lo0']*data.nvisit:(1 + data.par_order['lo0'])*data.nvisit]
 	print "size\t", m.params[data.par_order['size']*data.nvisit:(1 + data.par_order['size'])*data.nvisit]
 	print "spot_T\t", m.params[data.par_order['spot_T']*data.nvisit:(1 + data.par_order['spot_T'])*data.nvisit]
 	print "p_T\t", m.params[data.par_order['p_T']*data.nvisit:(1 + data.par_order['p_T'])*data.nvisit]
 	print "sph0\t", m.params[data.par_order['sph0']*data.nvisit:(1 + data.par_order['sph0'])*data.nvisit]
 	print "sph1\t", m.params[data.par_order['sph1']*data.nvisit:(1 + data.par_order['sph1'])*data.nvisit]
 	print "sph2\t", m.params[data.par_order['sph2']*data.nvisit:(1 + data.par_order['sph2'])*data.nvisit]
 	print "sph3\t", m.params[data.par_order['sph3']*data.nvisit:(1 + data.par_order['sph3'])*data.nvisit]

class Model:
	"""
	doc
	"""
	def __init__(self, params, data, flags):
		p = FormatParams(params, data)
		self.lc = np.zeros(len(data.time))				#full model light curve (with systematics)
		self.transit_model = np.zeros(len(data.time))			#transit model (relative flux; no systematics)
		self.eclipse_model = np.zeros(len(data.time))			#eclipse model (relative flux; no systematics)
		self.phase_model = np.zeros(len(data.time))			#phase curve model (relative flux; no systematics; includes eclipse)
		self.phase_model_no_eclipse = np.zeros(len(data.time))		#phase variation model (relative flux; no systematics; no eclipse))
		self.data_corr = np.zeros(len(data.time))			#data with the odd/even effect and orbit-long ramps removed
		self.phase = np.zeros(len(data.time))				#orbital phase	(defined in model because it depends on the period and ephemeris)
		self.sys = np.zeros(len(data.time))				#systematics showing visit-long trends only
		self.all_sys = np.zeros(len(data.time))				#all systematics 
		self.data_normalized = np.zeros(len(data.time))
		self.bestfit = np.zeros(len(data.time))
		self.bestfit_no_eclipse = np.zeros(len(data.time))

		for i in range(data.nvisit):
			ind = data.vis_num == i
			self.phase[ind] = (data.time[ind] - p.t0[i])/p.per[i] - np.floor((data.time[ind] - p.t0[i])/p.per[i])
			if data.lc_type == "eclipse": 
				self.lc[ind] = get_elc(data.time[ind], p, data, i) 	
				#print "need to add dilution here"
			elif data.lc_type == "transit": 
				self.lc[ind] = get_tlc(data.time[ind], p, data, i) + data.dilution	#LK addition
			elif data.lc_type == "phase_curve":
				self.eclipse_model[ind] = get_elc(data.time[ind], p, data, i) 	
				self.transit_model[ind] = get_tlc(data.time[ind], p, data, i) 
				self.phase_model[ind] = get_phaselc(data.time[ind], p, data, i)
				self.lc[ind] = self.transit_model[ind] + (self.eclipse_model[ind] - 1.)*self.phase_model[ind]
				self.bestfit[ind] = self.transit_model[ind] + (self.eclipse_model[ind] - 1.)*self.phase_model[ind]
				self.bestfit_no_eclipse[ind] = self.transit_model[ind] + (np.max(self.eclipse_model) - 1.)*self.phase_model[ind]
				#print "LK: need to correct for dilution here"
			elif data.lc_type == "physical_model":
 				self.transit_model[ind] = get_tlc(data.time[ind], p, data, i) 
 				#self.phase_model[ind] = get_spidermanzhang(data.time[ind], p, data, i, data.stellar_grid)
 				#self.phase_model[ind] = get_spidermanspherical(data.time[ind], p, data, i, data.stellar_grid)
 				self.phase_model[ind] = get_spidermanhotspot(data.time[ind], p, data, i, data.stellar_grid)
				print "choosing spiderman model by hand here2"
 				#self.lc[ind] = (self.transit_model[ind]-1.0) + self.phase_model[ind]
 				self.lc[ind] = (self.transit_model[ind]-1.0) + (self.phase_model[ind]-1.)/(1. + data.dilution) + 1.
				self.bestfit[ind] = self.lc[ind]
 			else: assert False, "Unknown option; supported light curve types are 'transit', 'eclipse', 'phase_curve' and 'physical_model'"

			if flags['divide-white'] == False:
				S = np.ones_like(self.transit_model[ind])+p.scale[i]*data.scan_direction[ind]
				D = np.ones_like(self.transit_model[ind])+p.r3[i]*data.t_delay[ind]
				#self.all_sys[ind] = data.flux[ind]/self.lc[ind]
				#self.all_sys[ind]  = (p.c[i]*S + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*(1.0-np.exp((-p.r1[i]*data.t_orb[ind]-p.r2[i]-D)))
				self.all_sys[ind]  = (p.c[i] + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*(1.0-np.exp((-p.r1[i]*data.t_orb[ind]-p.r2[i]-D)))
				self.lc[ind] *= (p.c[i]*S + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*(1.0-np.exp((-p.r1[i]*data.t_orb[ind]-p.r2[i]-D)))
				self.data_corr[ind] = data.flux[ind]/((p.c[i]*S + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*(1.0-np.exp(-p.r1[i]*data.t_orb[ind]-p.r2[i]-D)))
				self.sys[ind] = p.c[i]+ p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2
				self.data_normalized[ind] = data.flux[ind]/(1.-np.exp((-p.r1[i]*data.t_orb[ind]-p.r2[i]-D))) - self.transit_model[ind]*p.c[i]*(S-1.)
			else:
				self.lc[ind] *= (p.c[i] + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*data.all_sys[ind]
				self.data_corr[ind] = data.flux[ind]/((p.c[i] + p.v[i]*data.t_vis[ind] + p.v2[i]*data.t_vis[ind]**2)*data.all_sys[ind])

		self.resid = data.flux - self.lc 
		self.chi2 = np.sum((self.resid/data.err)**2)		
		self.chi2red = self.chi2/data.dof
		self.rms = 1.0e6*np.sqrt(np.mean((self.resid/data.flux)**2))
		self.rms_predicted = 1.0e6*np.sqrt(np.mean((data.err/data.flux)**2))
		self.data = data
		self.ln_likelihood = -0.5*(np.sum((self.resid/data.err)**2 + np.log(2.0*np.pi*(data.err)**2)))
		self.bic = -2.*self.ln_likelihood + data.nfree_param*np.log(data.npoints)

def plot(params, data, flags, obs_par, plot_sys=False):
	p = FormatParams(params, data)
	m = Model(params,data, flags)
	sns.set_palette("muted")
	palette = sns.color_palette("muted", m.data.nvisit)
	plt.clf()
	if plot_sys==False:
		for i in range(m.data.nvisit):
			ind = m.data.vis_num == i
			plt.subplot(211)
			if data.lc_type == "transit":
				temp_phase = np.copy(m.phase[ind])
				temp_phase[temp_phase>0.5] -= 1.0
			 	plt.plot(temp_phase, m.data_corr[ind], marker='o', linestyle="None", markersize=5)
			else: plt.plot(m.phase[ind], m.data_corr[ind], marker='o', linestyle="None", markersize=5)
			plt.subplot(212)
			if data.lc_type == "transit": plt.plot(temp_phase, 1.0e6*m.resid[ind]/m.data.flux[ind], marker='o', linestyle="None", markersize=5)
			else: plt.plot(m.phase[ind], 1.0e6*m.resid[ind]/m.data.flux[ind], marker='o', linestyle="None", markersize=5)
		plt.subplot(211)
		phase_hr = np.linspace(m.phase.min()-0.05, m.phase.max()+0.05, 1000)
		t_hr = phase_hr*p.per[0]+p.t0[0]

		if data.lc_type == "transit":
			transit_model_hr = get_tlc(t_hr, p, data, 0) + data.dilution
			lc_hr = transit_model_hr
		elif data.lc_type == "eclipse":
			eclipse_model_hr = get_elc(t_hr, p, data, 0) 
			lc_hr = eclipse_model_hr
	#		print "need to add dilution here"
		elif data.lc_type == "phase_curve":
			transit_model_hr = get_tlc(t_hr, p, data, 0) 
			eclipse_model_hr = get_elc(t_hr, p, data, 0) 
			phase_model_hr = get_phaselc(t_hr, p, data, 0)
			lc_hr = transit_model_hr + (eclipse_model_hr-1.)*phase_model_hr
			print "need to add dilution here"
		elif data.lc_type == "physical_model":
 			transit_model_hr = get_tlc(t_hr, p, data, 0) 
 			#physical_model_hr = get_spidermanzhang(t_hr, p, data, 0, data.stellar_grid)
 			#physical_model_hr = get_spidermanspherical(t_hr, p, data, 0, data.stellar_grid)
 			physical_model_hr = get_spidermanhotspot(t_hr, p, data, 0, data.stellar_grid)
			print "choosing spiderman model by hand here1"
 			#lc_hr = (transit_model_hr-1.0) + physical_model_hr
 			lc_hr = (np.array(transit_model_hr)-1.0) + (np.array(physical_model_hr) - 1.)/(1.+data.dilution) + 1.
 		else: assert False, "Unknown option; supported light curve types are 'transit', 'eclipse', 'phase_curve' and 'physical_model'"

		if data.lc_type == "transit": 
			phase_hr[phase_hr > 0.5] -= 1.0
			ind = np.argsort(phase_hr)
			phase_hr = phase_hr[ind]
			lc_hr = lc_hr[ind]
			
		plt.plot(phase_hr, lc_hr, color='0.2', zorder=1)
		plt.axvline(0.57)	#LK temp
		plt.axvline(0.72)	#LK temp
		plt.ylabel("Relative flux")
		ax = plt.gca()
		ax.text(0,1,'obs, exp rms (ppm); chi:\n '+'{0:d}'.format(int(m.rms))+", "+'{0:d}'.format(int(m.rms_predicted))+", "+'{0:0.2f}'.format(m.chi2red),verticalalignment='top',horizontalalignment='left',transform=ax.transAxes, fontsize=10)
		if flags['fit-white']: plt.title("Fit to white light curve")
		else: plt.title("Fit to {0:0.2f}".format(m.data.wavelength)+" micron channel")
		delta = 3.0*m.rms/1.0e6
		plt.ylim((lc_hr.min() - delta , lc_hr.max()+delta))
		if data.lc_type == "phase_curve": plt.ylim((1.0-delta, lc_hr.max()+delta))
		if data.lc_type == "physical_model": plt.ylim((1.0-delta, lc_hr.max()+delta))	
		plt.xlim((phase_hr.min(), phase_hr.max()))
		if obs_par['lc_type'] == 'transit': plt.xlim(-0.2, 0.2)

		plt.subplot(212)
		plt.axhline(0, zorder=1, color='0.2', linestyle='dashed')
		plt.ylabel("Residuals (ppm)")
		plt.xlabel("Orbital phase")
		plt.xlim((phase_hr.min(), phase_hr.max()))
		if obs_par['lc_type'] == 'transit': plt.xlim(-0.2, 0.2)

		if flags['output']: 
			figname =  "fit" + "{0:0.2f}".format(data.wavelength) + ".png"
			plt.savefig(figname)
		plt.show()

	elif plot_sys == True:
		#offset = np.mean(m.sys)*0.001
		print "FIXME: this option isn't implemented for phase curves"
		offset = 0.
		for i in range(m.data.nvisit):
			ind = m.data.vis_num == i
			plt.subplot(211)
			plt.plot(m.phase[ind], m.data_normalized[ind]+offset*i, marker='o', linestyle="None", markersize=5)
			plt.plot(m.phase[ind], m.sys[ind]*m.transit_model[ind]+offset*i, color=palette[i], zorder=1)
			plt.plot(m.phase[ind], m.sys[ind]+offset*i, zorder=1, linestyle='dashed')
			plt.ylabel("Corrected flux")
			plt.title("Flux corrected for orbit-long trends and scan direction")
			plt.subplot(212)
			plt.plot(m.phase[ind], m.resid[ind], marker='o', linestyle="None", markersize=5, color=palette[i])
			plt.xlabel("Orbital phase")
			plt.ylabel("Residuals (e-)")
		plt.subplot(211)
		plt.plot(m.phase[0], m.sys[0]*m.transit_model[0], color='0.5', zorder=1, label="with transit")
		plt.plot(m.phase[0], m.sys[0], color='0.5', zorder=1, linestyle='dashed', label="transit removed")
		plt.legend()
		plt.show()
			
def usage():
	cmd = sys.argv[0]
	sys.stderr.write('Usage: python %s OPTION\n\n' % os.path.basename(cmd))
	sys.stderr.write(
		'Allowed OPTION flags are:\n'
		'  --show-plot      		displays fitted light curve plots\n'
		'  --run-mcmc      		runs MCMC starting from least-squares best fit parameters\n'
		'  --run-pb         		runs prayer bead analysis\n'
		'  --plot-raw-data		plots raw light curve separated by visit\n'   
		'  --plot-sys			plots light curve fit with visit-long systematics included\n'   
		'  --path PATHNAME		specifies PATHNAME to light curves to be fit (default = ./spec_lc/*)\n'
		'  --fit-white FILENAME		fits the white light curve stored in FILENAME (default = "lc_white.txt"\n'
		'  -v               		prints fit diagnostic information\n'
		'  -o               		saves fit output to file\n'
		'  -h               		lists instructions for usage\n'
		'\n')
	sys.exit(1)

def make_dict(table):
	return {x['parameter']: x['value'] for x in table}

def plot_data(data):
	#palette = sns.color_palette("husl", data.nvisit)
	for i in range(data.nvisit): 	
		ind = data.vis_num==i
		plt.subplot((data.nvisit)*100+10+i+1)
		plt.plot(data.t_vis[ind]*24., data.flux[ind], marker='o', markersize=4.5, linestyle="none", label = "Visit {0}".format(i))
		plt.xlim(((data.t_vis.min()-0.02)*24., (data.t_vis.max()+0.05)*24.))
		plt.ylim((0.998*data.flux.min(), 1.002*data.flux.max()))
		#ind = data.orb_num == 3
		#plt.plot(data.t_vis[ind]*24., data.flux[ind], 'xk')
		plt.legend()
	plt.xlabel("Time after visit start (hours)")
	plt.ylabel("Flux (e-)")
	plt.tight_layout()
	plt.show()	


def weighted_mean(data, err):            #calculates the weighted mean for data points data with std devs. err
	ind = err != 0.0
        weights = 1.0/err[ind]**2
        mu = np.sum(data[ind]*weights)/np.sum(weights)
        var = 1.0/np.sum(weights)
        return [mu, np.sqrt(var)]                

def get_tlc(t, p, data, v_num):
	bat_params = batman.TransitParams()
	bat_params.t0 = p.t0[v_num]
	bat_params.t_secondary = p.t_secondary[v_num]
	bat_params.per = p.per[v_num]
	bat_params.rp = p.rp[v_num]
	bat_params.a = p.a[v_num]
	bat_params.inc = p.inc[v_num]
	bat_params.ecc = p.ecc[v_num]
	bat_params.w = p.w[v_num]
	bat_params.u = [p.u1[v_num], p.u2[v_num]]
	bat_params.limb_dark = "quadratic"		#FIXME - specify this value in one of the config files
	
	m = batman.TransitModel(bat_params, t, supersample_factor=3, exp_time = data.exp_time/24./60./60.)
	return m.light_curve(bat_params)

def get_elc(t, p, data, v_num):
	bat_params = batman.TransitParams()
	bat_params.t0 = p.t0[v_num]
	bat_params.t_secondary = p.t_secondary[v_num]
	bat_params.per = p.per[v_num]
	bat_params.rp = p.rp[v_num]
	bat_params.a = p.a[v_num]
	bat_params.inc = p.inc[v_num]
	bat_params.ecc = p.ecc[v_num]
	bat_params.w = p.w[v_num]
	bat_params.u = [p.u1[v_num], p.u2[v_num]]
	bat_params.fp = p.fp[v_num] 
	bat_params.limb_dark = "quadratic"		#FIXME - specify this value in one of the config files
	
	m = batman.TransitModel(bat_params, t, supersample_factor=3, exp_time = data.exp_time/24./60./60., transittype="secondary")
	return m.light_curve(bat_params)

def get_phaselc(t, p, data, v_num):
	return 1.+p.amp1[v_num]*np.cos(2.*np.pi*(t-p.theta1[v_num])/p.per[v_num]) + p.amp2[v_num]*np.cos(4.*np.pi*(t-p.theta2[v_num])/p.per[v_num])


def proj_area(phi, inc):
        R0 = 1.07e8             #planet radius in m
        p = 0.00117             # planet-to-star mass ratio
        r = 2.9695e9            #orbital separation in m
	kn = 0.653		#from table B5 of Leconte et al. 2011
        n = 1.005
        qn = kn*(1.- n/5.)      #n 


        alpha1 = 2.5*qn*R0**3/(p*r**3)
        alpha2 = -1.25*qn*R0**3/(p*r**3)
        alpha3 = alpha2

        a1 = R0*(1+alpha1)
        a2 = R0*(1+alpha2)
        a3 = R0*(1+alpha3)

        #plt.plot(phi, np.sqrt(a3**2*np.sin(inc)**2*(a1**2*np.sin(phi)**2+a2**2*np.cos(phi)**2)+ a1**2*a2**2*np.cos(inc)**2)/(a2*a3), '.k')
	#plt.show()
        return  np.sqrt(a3**2*np.sin(inc)**2*(a1**2*np.sin(phi)**2+a2**2*np.cos(phi)**2)+ a1**2*a2**2*np.cos(inc)**2)/(a2*a3)

def get_spidermanzhang(t, p, data, v_num, stellar_grid):
        web_p = spiderman.ModelParams(brightness_model =  'zhang') 
        
        web_p.n_layers = 5
        web_p.t0 = p.t0[v_num] 
        web_p.per = p.per[v_num]
        web_p.a_abs = p.a_abs[v_num]
        web_p.inc = p.inc[v_num]
        web_p.ecc = p.ecc[v_num]
        web_p.w = p.w[v_num] 
        #web_p.rp = p.rp[v_num] 
        web_p.rp = 0.1127
        web_p.a = p.a[v_num] 
        web_p.p_u1 = p.p_u1[v_num] 
        web_p.p_u2 = p.p_u2[v_num]
        web_p.T_s = p.T_s[v_num]
        web_p.l1 = data.l1
        web_p.l2 = data.l2
	web_p.xi = p.xi[v_num]
	web_p.T_n = p.T_n[v_num]
	web_p.delta_T = p.delta_T[v_num]
	web_p.delta_T_cloud = p.delta_T_cloud[v_num]
	web_p.thermal = True

        lc = spiderman.web.lightcurve(t, web_p, stellar_grid = stellar_grid)

	phs = (t - p.t0[v_num])/p.per[v_num]
	phs -= np.round(phs)
	rprs2 = proj_area(phs*2.*np.pi, p.inc[v_num]*np.pi/180.)

	#print np.max(rprs2), np.min(rprs2)
	return (np.array(lc) - 1.0)*rprs2 + 1.

def get_spidermanspherical(t, p, data, v_num, stellar_grid):
        web_p = spiderman.ModelParams(brightness_model =  'spherical', thermal = True) 
        
        web_p.n_layers = 5
        web_p.t0 = p.t0[v_num] 
        web_p.per = p.per[v_num]
        web_p.a_abs = p.a_abs[v_num]
        web_p.inc = p.inc[v_num]
        web_p.ecc = p.ecc[v_num]
        web_p.w = p.w[v_num] 
        #web_p.rp = p.rp[v_num] 
        web_p.rp = 0.1127
        web_p.a = p.a[v_num] 
        web_p.p_u1 = p.p_u1[v_num] 
        web_p.p_u2 = p.p_u2[v_num]
        web_p.T_s = p.T_s[v_num]
        web_p.l1 = data.l1
        web_p.l2 = data.l2
	web_p.degree = 2
	web_p.la0 = p.la0[v_num]
	web_p.lo0 = p.lo0[v_num]
	sph0 = p.sph0[v_num]
	sph1 = p.sph1[v_num]
	sph2 = p.sph2[v_num]
	sph3 = p.sph3[v_num]
	web_p.sph = [sph0, sph1, sph2, sph3]

        lc = spiderman.web.lightcurve(t, web_p, stellar_grid = stellar_grid)

	phs = (t - p.t0[v_num])/p.per[v_num]
	phs -= np.round(phs)
	rprs2 = proj_area(phs*2.*np.pi, p.inc[v_num]*np.pi/180.)

	return (np.array(lc) - 1.0)*rprs2 + 1.

def get_spidermanhotspot(t, p, data, v_num, stellar_grid): 
        web_p = spiderman.ModelParams(brightness_model =  'hotspot_t') 
        
        web_p.n_layers = 5
        web_p.t0 = p.t0[v_num] 
        web_p.per = p.per[v_num]
        web_p.a_abs = p.a_abs[v_num]
        web_p.inc = p.inc[v_num]
        web_p.ecc = p.ecc[v_num]
        web_p.w = p.w[v_num] 
        #web_p.rp = p.rp[v_num] 
        web_p.rp = 0.1127
        web_p.a = p.a[v_num] 
        web_p.p_u1 = p.p_u1[v_num] 
        web_p.p_u2 = p.p_u2[v_num]
        web_p.T_s = p.T_s[v_num]
        web_p.l1 = data.l1
        web_p.l2 = data.l2
	web_p.la0 = p.la0[v_num]
	web_p.lo0 = p.lo0[v_num]
	web_p.size = p.size[v_num]
	web_p.spot_T = p.spot_T[v_num]
	web_p.p_T = p.p_T[v_num]
	web_p.thermal = True
		
        lc = spiderman.web.lightcurve(t, web_p, stellar_grid = stellar_grid)

	phs = (t - p.t0[v_num])/p.per[v_num]
	phs -= np.round(phs)
	rprs2 = proj_area(phs*2.*np.pi, p.inc[v_num]*np.pi/180.)

	return (np.array(lc) - 1.0)*rprs2 + 1.

def residuals(params, data, flags, fjac=None):					
	return [0, Model(params, data, flags).resid/data.err]

def least_sq_fit(file_name, obs_par, fit_par, data, flags):
	nvisit = int(obs_par['nvisit'])
	npar = len(fit_par)*nvisit

	#initializes least squares fit parameters
	parinfo = [{'value':0, 'fixed':0, 'limited':[0,0,], 'limits':[0.0,0.0], 'step':0.0} for j in range(npar)]
	params_s = []

	for i in range(npar/nvisit):								#loops through params
		for j in range(nvisit):								#loops through visits
			parinfo[i*nvisit+j]['value'] = fit_par['value'][i]			#sets initial guess value	
			parinfo[i*nvisit+j]['step'] = 0.01*np.abs(fit_par['value'][i])		#sets parameter step size
			parinfo[i*nvisit+j]['fixed'] = fit_par['fixed'][i].lower() == "true"	#sets whether parameter varies
			if j>0 and fit_par['tied'][i].lower() == "true":
				parinfo[i*nvisit+j]['tied'] = 'p[{0}]'.format(nvisit*i)		#ties parameters to first visit value
			if fit_par['lo_lim'][i].lower() == "true": 				#puts lower limit on parameter
				parinfo[i*nvisit+j]['limited'][0] = True
				parinfo[i*nvisit+j]['limits'][0] = fit_par['lo_val'][i]
			if fit_par['hi_lim'][i].lower() == "true": 				#puts upper limit on parameter
				parinfo[i*nvisit+j]['limited'][1] = True
				parinfo[i*nvisit+j]['limits'][1] = fit_par['hi_val'][i]
			if (i == 13):					
				parinfo[i*nvisit+2]['value'] = 0.
				parinfo[i*nvisit+2]['fixed'] = True
				parinfo[i*nvisit+3]['value'] = 0.
				parinfo[i*nvisit+3]['fixed'] = True
			params_s.append(fit_par['value'][i])
	print "fixing v2 at 0 for zhao data"
	#LK

	params_s = np.array(params_s)

	if flags['plot-raw-data']: plot_data(data)
	fa = {'data':data, 'flags':flags}

	if flags['divide-white']:
		sys_vector = np.genfromtxt("white_systematics.txt")
		data.all_sys = sys_vector

	m = mpfit.mpfit(residuals, params_s, functkw=fa, parinfo = parinfo, quiet=True) 
	model = Model(m.params, data, flags)

	#LK test:	give outliers zero weight in fit
	ind = model.resid/data.err>4.
	print "Number of outliers", sum(ind)
	data.err[ind] = 100000000.
	m = mpfit.mpfit(residuals, params_s, functkw=fa, parinfo = parinfo, quiet=True) 
	model = Model(m.params, data, flags)
	
	if flags['output']: 
		f = open(flags['out-name'], "a")
		print>>f, "{0:0.3f}".format(data.wavelength), "{0:0.6f}".format(m.params[data.par_order['fp']*nvisit]), "{0:0.6f}".format(m.perror[data.par_order['fp']*nvisit]), "{0:0.6f}".format(m.params[data.par_order['rp']*nvisit]), "{0:0.6f}".format(m.perror[data.par_order['rp']*nvisit]), "{0:0.2f}".format(Model(m.params, data, flags).chi2red), "{0:0.3f}".format(m.params[data.par_order['u1']*nvisit]), "{0:0.3f}".format(m.perror[data.par_order['u1']*nvisit])
		f.close()

		edepth = m.params[data.par_order['fp']*nvisit]
		edeptherr = m.perror[data.par_order['fp']*nvisit]
		pickle.dump([data.wavelength, model.phase, model.data_corr, data.err/data.flux, model.resid/data.flux, edepth, edeptherr, model.lc, data.flux, data.vis_num], open("fcorr"+"{0:0.2f}".format(data.wavelength)+".p", "wb"))
		pickle.dump([data, model, m.params], open("bestfit_"+"{0:0.2f}".format(data.wavelength)+".pic", "wb"),-1)

	if flags['verbose']: 
		model = Model(m.params, data, flags)
		print "wavelength:", "{0:0.3f}".format(data.wavelength), "\tchi2_red", "{0:0.2f}".format(model.chi2red), "\trms", "{0:0f}".format(model.rms), "\tBIC", "{0:0.2f}".format(model.bic), "\tnfree_param", data.nfree_param
		PrintParams(m, data)

	if flags['show-plot']: plot(m.params, data, flags, obs_par, plot_sys=flags['plot-sys'])
	return m.params

def format_params_for_mcmc(params, obs_par, fit_par):	#FIXME: make sure this works for cases when nvisit>1
	nvisit = int(obs_par['nvisit'])				
	theta = []

	for i in range(len(fit_par)):
		if fit_par['fixed'][i].lower() == "false":
			if fit_par['tied'][i].lower() == "true": theta.append(params[i*nvisit])
			else: 
			#	for j in range(nvisit): theta.append(params[i*nvisit+j])
				if(i != 13): 
					for j in range(nvisit): theta.append(params[i*nvisit+j])
				else: 
					for j in range(nvisit -2): theta.append(params[i*nvisit+j])
				#fix v2 for zhao eclipses by hand
				#LK 
	return np.array(theta)

def format_prior_for_mcmc(obs_par, fit_par):
	nvisit = int(obs_par['nvisit'])				
	prior = []

	for i in range(len(fit_par)):
		if fit_par['fixed'][i].lower() == "false":
			if fit_par['tied'][i].lower() == "true": prior.append([fit_par['prior'][i], float(fit_par['p1'][i]), float(fit_par['p2'][i])])
			else: 
				#for j in range(nvisit): prior.append([fit_par['prior'][i], float(fit_par['p1'][i]), float(fit_par['p2'][i])])
				if(i != 13): 
					for j in range(nvisit): prior.append([fit_par['prior'][i], float(fit_par['p1'][i]), float(fit_par['p2'][i])])
				else: 
					for j in range(nvisit-2): prior.append([fit_par['prior'][i], float(fit_par['p1'][i]), float(fit_par['p2'][i])])
				#fix v2 for zhao eclipses by hand
				#LK
	return prior


def mcmc_output(samples, params, obs_par, fit_par, data, chain):	#FIXME: make sure this works for cases when nvisit>1
	nvisit = int(obs_par['nvisit'])				
	labels = []

	for i in range(len(fit_par)):
		if fit_par['fixed'][i].lower() == "false":
			if fit_par['tied'][i].lower() == "true": labels.append(fit_par['parameter'][i])
			else: 
				for j in range(nvisit): labels.append(fit_par['parameter'][i]+str(j))
	
	chainname = "chain_" + pythontime.strftime("%Y_%m_%d_%H:%M:%S") + str(data.wavelength)
	flatchainname = "flatchain_" + pythontime.strftime("%Y_%m_%d_%H:%M:%S") + str(data.wavelength)
	figname = "pairs_" + pythontime.strftime("%Y_%m_%d_%H:%M:%S") + str(data.wavelength) + ".png"

	np.save(chainname, chain)
	np.save(flatchainname, samples)
	#fig = corner.corner(samples, labels=labels, show_titles=True)
	#fig.savefig(figname)


def format_params_for_Model(theta, params, obs_par, fit_par):
	nvisit = int(obs_par['nvisit'])
	params_updated = []
	iter = 0									#this should be a more informative name FIXME
	for i in range(len(fit_par)):
		if fit_par['fixed'][i].lower() == "true": 
			for j in range(nvisit): 
				params_updated.append(params[i*nvisit+j])
		else:
			if fit_par['tied'][i].lower() == "true": 
				for j in range(nvisit): params_updated.append(theta[iter])
				iter += 1
			else: 
				#for j in range(nvisit): 		
				#	params_updated.append(theta[iter])
				#	iter += 1
				if(i != 13): 
					for j in range(nvisit): 
						params_updated.append(theta[iter])
						iter += 1
				else: 
					for j in range(nvisit-2): 
						params_updated.append(theta[iter])
						iter += 1
					params_updated.append(0.) 	#adds a zero for v2[2]
					params_updated.append(0.) 	#adds a zero for v2[3]
				#LK
	return np.array(params_updated)
		

def mcmc_fit(file_name, obs_par, fit_par, flags):
	data = LightCurveData(file_name, obs_par, fit_par, flags)
	if flags['divide-white']:
		sys_vector = np.genfromtxt("white_systematics.txt")
		data.all_sys = sys_vector


	params = least_sq_fit(file_name, obs_par, fit_par, data, flags)		#starting guess
	theta = format_params_for_mcmc(params, obs_par, fit_par)	

	ndim, nwalkers = len(theta), 100					#FIXME set nwalkers is a config file
	sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (params, data, obs_par, fit_par, flags), threads=14)

	pos = [theta + 1e-5*np.random.randn(ndim) for i in range(nwalkers)]

	#sampler.run_mcmc(pos,200)
	#sampler.run_mcmc(pos,5000)

	nsteps = 6000
	#nsteps = 3000
	for i, result in enumerate(sampler.sample(pos, iterations=nsteps)):
		#if i%100==0: print "% complete:", float(i)/float(nsteps) 
		print "step", i
		#if i > 5100:
		#	GR

	samples = sampler.chain[:, 4000:, :].reshape((-1, ndim))
	#samples_firsthalf = sampler.chain[:, 4000:5000, :].reshape((-1, ndim))
	#samples_secondhalf = sampler.chain[:, 5000:6000, :].reshape((-1, ndim))
	

	#samples = sampler.chain[:, 600:, :].reshape((-1, ndim))

	mcmc_output(samples, params, obs_par, fit_par, data, sampler.chain)

	medians = []
	errors = []

	for i in range(len(theta)):
		q = quantile(samples[:, i], [0.16, 0.5, 0.84])
		medians.append(q[1])
		errors.append(q[2] - q[1])
	return data.wavelength, medians[0], errors[0], samples


def lnprior(theta, params, data, obs_par, fit_par, flags):
	lnprior_prob = 0.
	n = len(data.prior)
	for i in range(n):
		if data.prior[i][0] == 'U': 
	#		print theta[i], data.prior[i][1], data.prior[i][2]
			if np.logical_or(theta[i] < data.prior[i][1], theta[i] > data.prior[i][2]): lnprior_prob += - np.inf
		if data.prior[i][0] == 'N': 
			lnprior_prob -= 0.5*(np.sum(((theta[i] - data.prior[i][1])/data.prior[i][2])**2 + np.log(2.0*np.pi*(data.prior[i][2])**2)))
		#if theta[19] > theta[17]: lnprior_prob += -np.inf
	return lnprior_prob
	

def lnprob(theta, params, data, obs_par, fit_par, flags):
	updated_params = format_params_for_Model(theta, params, obs_par, fit_par)
	lp = lnprior(theta, params, data, obs_par, fit_par, flags)
	if lp == -np.inf:
                return -np.inf
        else:
                m = Model(updated_params, data, flags)
                return m.ln_likelihood + lp


def main():
	#parses command line input
	try: opts, args = getopt.getopt(sys.argv[1:], "hov", ["help", "show-plot", "run-mcmc", "plot-raw-data", "plot-sys", "path=", "fit-white=", "divide-white"]) 
	except getopt.GetoptError: usage()

	#defaults for command line flags
	verbose, output, show_plot, run_mcmc, run_lsq, plot_raw_data, plot_sys, path, fit_white, divide_white = False, False, False, False, True, False, False, "spec_lc", False, False

	for o, a in opts:
		if o in ("-h", "--help"): usage()
		elif o == "-o": output = True
		elif o == "-v": verbose = True
		elif o == "--show-plot": show_plot = True
		elif o == "--run-mcmc": run_mcmc, run_lsq = True, False
		elif o == "--run-lsq": run_lsq = True
		elif o == "--plot-raw-data": plot_raw_data = True
		elif o == "--plot-sys": plot_sys = True
		elif o == "--path": path = a
		elif o == "--fit-white": fit_white, white_file = True, a
		elif o == "--divide-white": divide_white = True
		else: assert False, "unhandled option"

	flags = {'verbose': verbose, 'show-plot': show_plot, 'plot-raw-data': plot_raw_data, 'plot-sys': plot_sys, 'output': output, 'out-name': 'none.txt', 'run-lsq': run_lsq, 'run-mcmc': run_mcmc, 'divide-white': divide_white, 'fit-white': fit_white}		#put these in dictionary right away!!

	#reads in observation and fit parameters
	obs_par = make_dict(ascii.read("config/obs_par.txt", Reader=ascii.CommentedHeader))
	fit_par = ascii.read("config/fit_par.txt", Reader=ascii.CommentedHeader)		

	files = glob.glob(os.path.join(path, "*"))		
	if fit_white: files = glob.glob(white_file)

	flags['out-name'] = "analysis/fit_" + pythontime.strftime("%Y_%m_%d_%H:%M") + ".txt"

	if not flags['run-mcmc']:
		for f in files:
			data = LightCurveData(f, obs_par, fit_par, flags)
			if flags['divide-white']:
				sys_vector = np.genfromtxt("white_systematics.txt")
				data.all_sys = sys_vector
			params = least_sq_fit(f, obs_par, fit_par, data, flags)
			m = Model(params, data, flags)
	
		#outfile = open("white_systematics.txt", "w")
		#for i in range(len(m.all_sys)): print>>outfile, m.all_sys[i]
		#outfile.close()


	print "check to make sure you are doing the right thing for priors on delta_T_cloud"
	if flags['run-mcmc']:
		for f in files:
			output = mcmc_fit(f, obs_par, fit_par, flags)
			temp = "mcmc_out_"+'{0:0.2f}'.format(output[0])
			np.save(temp, output[3])
			print output[0], output[1], output[2]
			outfile = open("mcmc_output.txt", "a")
			print>>outfile, output[0], output[1], output[2]
			outfile.close()
				

if __name__ == '__main__':
	main()
