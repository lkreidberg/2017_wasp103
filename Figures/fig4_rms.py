import numpy as np
import matplotlib.pyplot as plt
from lc_fit_nospiderman import Model, LightCurveData
import pickle
import glob, os
import matplotlib.gridspec as gridspec
from matplotlib import rc
import seaborn as sns

sns.set_context("talk", font_scale = 1.0)
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

# COMPUTE ROOT-MEAN-SQUARE AND STANDARD ERROR OF DATA FOR VARIOUS BIN SIZES
def computeRMS(data, maxnbins=None, binstep=1, isrmserr=False):
    #data    = fit.normresiduals
    #maxnbin = maximum # of bins
    #binstep = Bin step size

    # bin data into multiple bin sizes
    npts    = data.size
    if maxnbins is None:
        maxnbins = npts/12.
    binsz   = np.arange(1, maxnbins+binstep, step=binstep)
    nbins   = np.zeros(binsz.size)
    rms     = np.zeros(binsz.size)
    rmserr  = np.zeros(binsz.size)
    for i in range(binsz.size):
        nbins[i] = int(np.floor(data.size/binsz[i]))
        bindata   = np.zeros(int(nbins[i]), dtype=float)
        # bin data
        # ADDED INTEGER CONVERSION, mh 01/21/12
        for j in range(int(nbins[i])):
            bindata[j] = data[int(j*binsz[i]):int((j+1)*binsz[i])].mean()
        # get rms
        rms[i]    = np.sqrt(np.mean(bindata**2))
        rmserr[i] = rms[i]/np.sqrt(2.*int(nbins[i]))
    # expected for white noise (WINN 2008, PONT 2006)
    stderr = (data.std()/np.sqrt(binsz))*np.sqrt(nbins/(nbins - 1.))
    if isrmserr == True:
        return rms, stderr, binsz, rmserr
    else:
        return rms, stderr, binsz


# Compute standard error
def computeStdErr(datastd, datasize, binsz):
    #datastd  = fit.normresiduals.std()
    #datasize = fit.normresiduals.size
    #binsz    = array of bins

    nbins   = np.zeros(binsz.size)
    for i in range(binsz.size):
        nbins[i] = int(np.floor(datasize/binsz[i]))
    stderr = (datastd/np.sqrt(binsz))*np.sqrt(nbins/(nbins - 1.))
    return stderr


path = './WFC3_best_fits/spec_fits/'
files = glob.glob(os.path.join(path, "*.pic"))

plt.figure(figsize=(5,6))

gs = gridspec.GridSpec(4, 3, hspace = 0.2, wspace=0.3) 
ncol = 3

for ii, f in enumerate(files):
    row, col = int(np.floor(ii/ncol)), ii%ncol
    ax = plt.subplot(gs[row, col])

    p = pickle.load(open(f, 'rb'))
    d, m, par = p[0], p[1], p[2] 
    
    ind = d.err < 9.0e7                     #indices of outliers
    resid = m.resid[ind]
    resid = resid[0:452]			#just phase curve data
 
    rms, stderr, binsz, rmserr = computeRMS(resid, maxnbins = None, binstep = 1, isrmserr = True)

    normfactor = stderr[0]
    #normfactor = m.rms_predicted

    plt.loglog(binsz, rms/normfactor, color='black', lw=1.5, label='Fit RMS')    # our noise
    plt.loglog(binsz, stderr/normfactor, color='blue', ls='-', lw=2, label='Std. Err.') # expected noise

    #plt.axhline(m.rms_predicted/normfactor, linestyle = 'dashed')

    plt.gca().set_xticks([1, 10])
    plt.gca().set_xticklabels([1, 10])
    plt.gca().set_yticks([1, 0.1])
    plt.gca().set_yticklabels([1, 0.1])

    ax.tick_params(axis='both', which='major', pad=2)

    plt.xlim(0, binsz[-1]*2)
    plt.ylim(stderr[-1]/normfactor/2., stderr[0]/normfactor*2.)
    plt.text(3.5, 1.2, str(d.wavelength) + " $\mu$m", fontsize = 12)

#    if row > 2: plt.xlabel("Bin Size") 
    if col == 0 and row == 1: 
	plt.ylabel("Normalized RMS", fontsize = 14) 
	ax.yaxis.set_label_coords(-0.3,0.)

#add Spitzer
#files = ['Ch1_best_fits/2017-10-19_17:04-zhang/bestfit.pic', 'Ch2_best_fits/2017-10-19_16:16-zhang/bestfit.pic']
files = ["Ch1_best_fits/2018-02-07_14:24-zhang/bestfit.pic", "Ch2_best_fits/2018-02-07_12:02-zhang/bestfit.pic"]
waves = [3.6, 4.5]
colors = ['orange', 'red']

for i, f in enumerate(files):
    ax = plt.subplot(gs[int(np.floor((10+i)/ncol)), (10+i)%ncol])

    p = pickle.load(open(f, "rb"))
    rms, stderr, binsz = p[17], p[18], p[19]

    normfactor = stderr[0]
    plt.loglog(binsz, rms/normfactor, color='black', lw=1.5, label='Fit RMS')    # our noise
    plt.loglog(binsz, stderr/normfactor, color=colors[i], ls='-', lw=2, label='Std. Err.') # expected noise
    plt.xlim(0, binsz[-1]*2)
    plt.ylim(stderr[-1]/normfactor/2., stderr[0]/normfactor*2.)
    plt.text(30, 1, str(waves[i]) + " $\mu$m", fontsize = 12)

    plt.gca().set_xticks([1, 10, 100])
    plt.gca().set_xticklabels([1, 10, 100])
    plt.gca().set_yticks([1, 0.1])
    plt.gca().set_yticklabels([1, 0.1])
    
    ax.tick_params(axis='both', which='major', pad=2)
    if i == 0: plt.xlabel("Points Per Bin", labelpad = 10, fontsize = 14)

plt.savefig("fig4.pdf")
plt.show()
