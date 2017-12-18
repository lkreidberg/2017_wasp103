import pickle
import numpy as np

import scipy.stats as st

def get_significance(alpha):
	z = st.norm.ppf(1.-alpha)
	return z

def quantile(x, q):
        return np.percentile(x, 100. * q)

p = pickle.load(open("Mike_models/WASP-103b_grid_DAYSIDE_output.pic", "rb"))
samples = p[9]


# metallicity
Z = samples[:,1]		#log metallicity
n = len(Z)*1.
alpha = sum(Z<0.)/n

print "significance that Z is > solar:", get_significance(alpha/2.)

###
#C/O
ctoo = samples[:,2]
n = len(ctoo)*1.0
ctoo = 10.**ctoo
ctoo_max = 0.9
alpha = sum(ctoo > ctoo_max)/n

print "significance that C/O <" , ctoo_max, ":", get_significance(alpha/2.)

print "C/O quantiles:", quantile(ctoo, np.array([0.16, 0.5, 0.84]))




