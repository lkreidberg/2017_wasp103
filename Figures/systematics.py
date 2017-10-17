import pickle
from astropy.io import ascii
import matplotlib.gridspec as gridspec
from pylab import *
from matplotlib import rc
import transit_lc
import seaborn as sns
from lc_fit import Model, LightCurveData
sns.set_context("talk", font_scale = 1.0)
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

gs = gridspec.GridSpec(3, 1, height_ratios =[1,1, 1],  hspace=0.00)

# plots HST white light phase curve
bestfits = ["./WFC3_best_fits/bestfit_zhang_allsys.pic"]

alpha = 0.7
ms = 5

plt.figure(figsize = (5, 7))

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

	scan = np.ones_like(phase)
	scan[d.scan_direction == 0] = 1. + par[d.par_order['scale']*d.nvisit]
	delta = 1. + par[d.par_order['scale']*d.nvisit]

	ind = d.vis_num == 0
	offset = 0.

	p0, = plt.plot(t[ind] - (t[ind])[0], delta*visit_sys[ind]+ offset, color = '0.7', label = 'systematics only', zorder=-10)
	p1, = ax.plot(t[ind] - (t[ind])[0], allsys[ind]*scan[ind] + offset, color = 'b', alpha = alpha)
	p2, = ax.plot(t[ind] - (t[ind])[0], d.flux[ind]*scan[ind] + offset, '.k', markersize = ms)

	ind = d.vis_num == 1
	plt.plot(t[ind] - (t[ind])[0], delta*visit_sys[ind]+ offset, color = '0.7', zorder = -10)
	p3, = ax.plot(t[ind] - (t[ind])[0], allsys[ind]*scan[ind] + offset, color = 'b', alpha = alpha)
	p4, = ax.plot(t[ind] - (t[ind])[0], d.flux[ind]*scan[ind] + offset, '.k', markersize=ms)


	ind = d.vis_num == 2
	offset = 1.44e7
	plt.plot(t[ind] - (t[ind])[0], delta*visit_sys[ind]+ offset, color = '0.7', zorder = -10)
	p5, = ax.plot(t[ind] - (t[ind])[0], allsys[ind]*scan[ind] + offset, color = 'b', alpha = alpha)
	p6, = ax.plot(t[ind] - (t[ind])[0], d.flux[ind]*scan[ind] + offset, '.k', markersize = ms)

	ind = d.vis_num == 3
	offset = 1.47e7
	plt.plot(t[ind] - (t[ind])[0], delta*visit_sys[ind]+ offset, color = '0.7', zorder = -10)
	p7, = ax.plot(t[ind] - (t[ind])[0], allsys[ind]*scan[ind] + offset, color = 'b', alpha = alpha) 
	p8, = ax.plot(t[ind] - (t[ind])[0], d.flux[ind]*scan[ind] + offset, '.k', markersize = ms)

	plt.text(1.0, 6.71e7, "HST/\nWFC3") 

plt.gca().set_xlabel([])
plt.gca().set_xticks([])
plt.ylim(6.7e7, 6.82e7)
plt.xlim(-0.1, 1.3)

plt.legend([(p1, p2), p0], ["Best fit", "systematics only"], frameon=True)

#plot Spitzer phase curves
colors = ['orange', 'red']
observation = ['Spitzer Ch. 1', 'Spitzer Ch. 2']

bestfits = ["Ch1_best_fits/2017-10-11_20:25-zhang/bestfit.pic", "Ch2_best_fits/2017-10-11_20:24-zhang/bestfit.pic"]

for ii, f in enumerate(bestfits):
	ax = plt.subplot(gs[ii+1, 0])
	
	p = pickle.load(open(f, 'rb'))
	data_corr = p[1] 
	err = p[2] 
	bestfit = p[3]
	phase = p[0] 

	abscissauc = p[10]
	binfluxuc = p[11]
	binstduc = p[12]
	bestfit = p[13]
	abscissa = p[15]
	sys = p[16]

	plt.plot(abscissa - abscissa[0], sys, color= '0.7', zorder=-10)
	plt.plot(abscissa - abscissa[0], bestfit, color = colors[ii], alpha = alpha)
	plt.errorbar(abscissauc - abscissa[0], binfluxuc, binstduc, fmt = '.k', markersize = ms)

	plt.text(0.8, 0.986*np.median(binfluxuc), observation[ii])

	plt.xlim(-0.1, 1.3)
	
	if ii == 1: plt.xlabel("Time since first exposure (days)")
	if ii == 0: 
		plt.ylabel("Photoelectrons (+ constant)")
		plt.gca().set_xticks([])
		plt.gca().set_xlabel([])

plt.tight_layout()
plt.savefig("systematics.pdf")
plt.show()
