import glob
import os
import numpy as np

def quantile(x, q):
        return np.percentile(x, 100. * q)

path = "/Users/lkreidberg/Desktop/Projects/Observations/HST/WASP103_HST_all/SINE_MODEL_MCMC"
files = glob.glob(os.path.join(path, "mcmc*npy"))

#equation for phase variation
#return 1.+p.amp1[v_num]*np.cos(2.*np.pi*(t-p.theta1[v_num])/p.per[v_num]) + p.amp2[v_num]*np.cos(4.*np.pi*(t-p.theta2[v_num])/p.per[v_num])
ind  = 17       #phase shift
offset = 61     #SINE_MODEL_MCMC
ind2  = 16      #amplitude
ind3 = 1        #fpfs

#ind  = 17      #phase shift
#offset = 2     #SINE_MODEL_MCMC_6bins

#ind = 15       #xi
#offset = 0     #PHYSICAL_MODEL_MCMC

for f in files:
        d = np.load(f)
        #print f[32+offset:36+offset], np.median(d[:,ind]), np.median(d[:,ind]) - quantile(d[:,ind], 0.16), quantile(d[:,ind], 0.84) - np.median(d[:,ind]), np.median(d[:,ind2]), np.median(d[:,ind2]) - quantile(d[:,ind2], 0.16), quantile(d[:,ind2], 0.84) - np.median(d[:,ind2]), np.median(d[:, ind3]), np.median(d[:,ind3]) - quantile(d[:,ind3], 0.16), quantile(d[:,ind3], 0.84) - np.median(d[:,ind3])
        print f[32+offset:36+offset], np.median(d[:,ind]), np.median(d[:,ind2]), np.median(d[:,ind3])
