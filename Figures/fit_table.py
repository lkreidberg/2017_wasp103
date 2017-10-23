import numpy as np
from lc_fit import Model, LightCurveData
import pickle

#zhang, hotspot, spherical

ch1_chi2 = [0.9, 0.9, 0.9]
chi1_aic = [7131.0, 7107.4, 7096.7]
chi1_bic = [7075.2, 7058.5, 7040.9]

ch2_chi2 = [1.07, 1.07277, 1.07215]
chi2_aic = [8345.95, 8355.5, 8331.6]
chi2_bic = [8401.61, 8384.2, 8387.3]


# bset fits for WFC3
wfc3_chi2 = []
wfc3_aic = []
wfc3_bic = []

models = ["zhang", "hotspot_t", "spherical"]
path = "WFC3_best_fits/"
for i, model in enumerate(models):

	p = pickle.load(open(path+"bestfit_"+model+".pic", "rb")) 

	d, m, par = p[0], p[1], p[2]            #stores data,  model, and best fit parameters into d, m, & par
	wfc3_bic.append(m.bic)
	wfc3_chi2.append(m.chi2)

print wfc3_bic
print wfc3_chi2
