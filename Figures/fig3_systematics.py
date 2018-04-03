import pickle
from astropy.io import ascii
import matplotlib.gridspec as gridspec
from pylab import *
from matplotlib import rc
import transit_lc
import seaborn as sns
from lc_fit import Model, LightCurveData
import os, glob
sns.set_context("talk", font_scale = 1.0)
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

gs = gridspec.GridSpec(6, 1, height_ratios =[1,1,1,1,1,1],  hspace=0.33)

# plots HST white light phase curve
#bestfits = ["./WFC3_best_fits/old_white_fits/bestfit_zhang_allsys.pic"]
bestfits = ["./WFC3_best_fits/bestfit_zhang_allsys.pic"]

alpha = 0.7
ms = 5

fig = plt.figure(figsize = (6.5, 10))

for ii, f in enumerate(bestfits):
	ax = plt.subplot(gs[ii, 0])

	p = pickle.load(open(f, 'rb'))
	d, m, par = p[0], p[1], p[2]		#stores data,  model, and best fit parameters into d, m, & par
	
	dilution = d.dilution + 1.

	ind = d.err < 9.0e7			#indices of outliers
	
	err = d.err[ind]/d.flux[ind]		#normalized per point uncertainty 
	phase = m.phase[ind]
	data_corr = m.data_corr[ind]
	allsys = m.lc[ind]	
	t = d.time[ind]
	visit_sys = m.all_sys[ind]

        data_corr /= 1e6
        allsys /= 1e6
        visit_sys /= 1e6
        d.flux /= 1e6

	print "WFC3 white: obs, exp rms (ppm)", np.std(m.resid[ind]/d.flux[ind])*1e6, np.sqrt(1./np.median(d.flux[ind]))*1e6

	scan = np.ones_like(phase)
	scan[d.scan_direction == 0] = 1. + par[d.par_order['scale']*d.nvisit]
	delta = 1. + par[d.par_order['scale']*d.nvisit]

	offset = 2457000.

	ind = d.vis_num == 0
	p0, = plt.plot(t[ind] - offset, delta*visit_sys[ind], color = '0.7', label = 'systematics', zorder=-10)
	p1, = ax.plot(t[ind] - offset, allsys[ind]*scan[ind] , color = 'k') 
	p2, = ax.plot(t[ind] - offset, d.flux[ind]*scan[ind] , '.b', markersize = ms, alpha = alpha)
        plt.legend([(p1, p2), p0], ["Best fit", "Systematics"], frameon=True, fontsize = 11)
	plt.text(0.03, 0.1, "HST/WFC3", fontsize = 12, transform=ax.transAxes, bbox=dict(boxstyle="round", ec="gray", fc="white"))
        plt.ylim(6.73e1, 6.77e1)

	ax = plt.subplot(gs[1, 0])
	ind = d.vis_num == 1
	plt.plot(t[ind] - offset, delta*visit_sys[ind], color = '0.7', zorder = -10)
	p3, = ax.plot(t[ind] - offset, allsys[ind]*scan[ind] , color = 'k') 
	p4, = ax.plot(t[ind] - offset, d.flux[ind]*scan[ind] , '.b', markersize=ms, alpha = alpha)
        plt.ylim(6.7e1, 6.74e1)
	plt.text(0.03, 0.1, "HST/WFC3", fontsize = 12, transform=ax.transAxes, bbox=dict(boxstyle="round", ec="gray", fc="white"))
	plt.gca().set_yticks(np.array([67.0, 67.2, 67.4]))


	ax = plt.subplot(gs[2, 0])
	ind = d.vis_num == 2
	plt.plot(t[ind] - offset, delta*visit_sys[ind], color = '0.7', zorder = -10)
	p5, = ax.plot(t[ind] - offset, allsys[ind]*scan[ind] , color = 'k') 
	p6, = ax.plot(t[ind] - offset, d.flux[ind]*scan[ind] , '.b', markersize = ms, alpha = alpha)
	plt.text(0.03, 0.1, "HST/WFC3", fontsize = 12, transform=ax.transAxes, bbox=dict(boxstyle="round", ec="gray", fc="white"))
	plt.ylabel("Photoelectrons ($\\times10^6$)")
        ax.yaxis.set_label_coords(-0.11, 0.1)

	ax = plt.subplot(gs[3, 0])
	ind = d.vis_num == 3
	plt.plot(t[ind] - offset, delta*visit_sys[ind], color = '0.7', zorder = -10)
	p7, = ax.plot(t[ind] - offset, allsys[ind]*scan[ind] , color = 'k') 
	p8, = ax.plot(t[ind] - offset, d.flux[ind]*scan[ind] , '.b', markersize = ms, alpha = alpha)

	plt.text(0.03, 0.1, "HST/WFC3", fontsize = 12, transform=ax.transAxes, bbox=dict(boxstyle="round", ec="gray", fc="white"))
	#plt.gca().set_yticks(np.array([6.72e7, 6.75e7, 6.78e7, 6.81e7]))
	#plt.gca().set_yticklabels(np.array(["6.72e7", "6.75e7", "6.78e7", "6.81e7"]))


#plt.gca().set_xlabel([])
#plt.gca().set_xticks([])
#plt.ylim(6.7e7, 6.82e7)
#plt.xlim(-0.1, 1.3)


#plot Spitzer phase curves
colors = ['orange', 'red']
observation = ['Spitzer Ch. 1', 'Spitzer Ch. 2']
depth = [4.5e-3, 5.7e-3]
ylo = [4.45, 2.39]
yhi = [4.56, 2.49]
fluxconv = [306.126, 266.648]
#calculated from Jonathan Fraine's code https://github.com/exowanderer/ExoplanetTSO/blob/master/ExoplanetTSO_Auxiliary.py
"""fluxConv  = testheader['FLUXCONV']
expTime   = testheader['EXPTIME']
gain      = testheader['GAIN']
fluxConversion = expTime*gain / fluxConv"""

toffsets = [2457170.7, 2457162.]

#bestfits = ["Ch1_best_fits/2017-10-11_20:25-zhang/bestfit.pic", "Ch2_best_fits/2017-10-11_20:24-zhang/bestfit.pic"]
bestfits = ["Ch1_best_fits/2018-02-07_14:24-zhang/bestfit.pic", "Ch2_best_fits/2018-02-07_12:02-zhang/bestfit.pic"]

for ii, f in enumerate(bestfits):
	ax = plt.subplot(gs[ii+4, 0])
	
	p = pickle.load(open(f, 'rb'))
	data_corr = p[1] 
	err = p[2] 
	bestfit = p[3]
	phase = p[0] 

        bjdtdb = p[8]

	abscissauc = p[10]
	binfluxuc = p[11]
	binstduc = p[12]
	bestfit = p[13]
	abscissa = p[15]
	sys = p[16]

        plt.plot(bjdtdb - offset, data_corr*fluxconv[ii]*np.mean(binfluxuc)/1e6, marker = '.', color = colors[ii], alpha = 0.2, markersize = 3, markeredgewidth = 0., linewidth = 0., zorder =-100000)
	plt.plot(abscissa - offset + toffsets[ii], sys*fluxconv[ii]/1e6, color= '0.7', zorder=-10)
	plt.plot(abscissa - offset + toffsets[ii], bestfit*fluxconv[ii]/1e6, color = 'k') 
	#plt.errorbar(abscissauc - abscissa[0], binfluxuc*fluxconv[ii]/1e6, binstduc*fluxconv[ii]/1e6, fmt = '.k', markersize = ms)
#	print (abscissa[1] - abscissa[0])*24.*60.	#time between exposures
	print "rms obs, exp (ppm)", 1.0e6*np.std((binfluxuc - bestfit)/binfluxuc), 1.e6/np.sqrt(np.median(binfluxuc*fluxconv[ii]))

	#plt.text(abscissa[0] - offset + toffsets[ii], 1.01*np.max(binfluxuc)*fluxconv[ii]/1e6, observation[ii], fontsize = 12, transform=ax.transAxes)
	plt.text(0.03, 0.1, observation[ii], fontsize = 12, transform=ax.transAxes, bbox=dict(boxstyle="round", ec="gray", fc="white"))


        if ii == 0: ax.set_xticks([170.8, 171.0, 171.2, 171.4, 171.6, 171.8])
	if ii == 1: plt.xlabel("Time (BJD_$\mathrm{TDB}$ - 2,457,000)", labelpad = 5)

plt.savefig("fig3.pdf")

#calculates rms for WFC3 spectroscopic best fits
path = "./WFC3_best_fits/spec_fits/"
files = glob.glob(os.path.join(path, "*"))	
for f in files: 
	p = pickle.load(open(f, 'rb'))
        d, m, par = p[0], p[1], p[2]            #stores data,  model, and best fit parameters into d, m, & par

        dilution = d.dilution + 1.

        ind = d.err < 9.0e7                     #indices of outliers

        print "WFC3 spec: obs, exp rms (ppm)", d.wavelength, m.rms, 1.0e6*np.sqrt(np.mean((d.err[ind]/d.flux[ind])**2)),  m.rms/(1.0e6*np.sqrt(np.mean((d.err[ind]/d.flux[ind])**2)))

