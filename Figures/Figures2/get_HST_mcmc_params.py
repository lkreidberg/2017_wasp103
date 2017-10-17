import glob
import os
import numpy as np
		
def quantile(x, q):
        return np.percentile(x, 100. * q)

#path = "../PHYSICAL_MODEL_MCMC/"
path = "../SINE_MODEL_MCMC"
#path = "../SINE_MODEL_MCMC_6bins"
files = glob.glob(os.path.join(path, "mcmc_out_*npy"))		
print files

ind  = 17	#phase shift
offset = -4 	#SINE_MODEL_MCMC
ind2  = 16	#amplitude
ind3 = 1	#fpfs

#ind  = 17	#phase shift
#offset = 2	#SINE_MODEL_MCMC_6bins

#ind = 15	#xi
#offset = 0	#PHYSICAL_MODEL_MCMC

for f in files:
	d = np.load(f)
	#print f[32+offset:36+offset], np.median(d[:,ind]), np.median(d[:,ind]) - quantile(d[:,ind], 0.16), quantile(d[:,ind], 0.84) - np.median(d[:,ind]), np.median(d[:,ind2]), np.median(d[:,ind2]) - quantile(d[:,ind2], 0.16), quantile(d[:,ind2], 0.84) - np.median(d[:,ind2])
	print f[32+offset:36+offset], np.median(d[:,ind2])*np.median(d[:,ind3])	#phase amplitude
	#need to do error propagation to get uncertainties
