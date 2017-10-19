import numpy as np
import matplotlib.pyplot as plt
from lc_fit_nospiderman import Model, LightCurveData
import pickle
import glob, os
import matplotlib.gridspec as gridspec

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

plt.figure(figsize=(5,7))

gs = gridspec.GridSpec(4, 3) 
ncol = 3

for ii, f in enumerate(files):
    ax = plt.subplot(gs[int(np.floor(ii/ncol)), ii%ncol])

    p = pickle.load(open(f, 'rb'))
    d, m, par = p[0], p[1], p[2] 
    
    ind = d.err < 9.0e7                     #indices of outliers
    resid = m.resid[ind]
    resid = resid[0:452]			#just phase curve data
 
    rms, stderr, binsz, rmserr = computeRMS(resid, maxnbins = None, binstep = 1, isrmserr = True)

    normfactor = stderr[0]
    plt.loglog(binsz, rms/normfactor, color='black', lw=1.5, label='Fit RMS')    # our noise
    plt.loglog(binsz, stderr/normfactor, color='red', ls='-', lw=2, label='Std. Err.') # expected noise
    plt.xlim(0, binsz[-1]*2)
    plt.ylim(stderr[-1]/normfactor/2., stderr[0]/normfactor*2.)
    plt.text(10, 1, str(d.wavelength) + " $\mu$m")

    #plt.xlabel("Bin Size", fontsize=14)
    #plt.ylabel("Normalized RMS", fontsize=14)


#add Spitzer



plt.tight_layout()
plt.show()
