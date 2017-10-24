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
	bic = m.chi2 + d.nfree_param*np.log(d.npoints)
	aic = m.chi2 + 2.*d.nfree_param
	print m.chi2red, bic, aic

#Temperatures:
#HST:
#model, Tmin, Tmax =  zhang 1942.83470024 3884.44683954
#model, Tmin, Tmax =  hotspot_t 0.0 2860.2228932
#model, Tmin, Tmax =  spherical 1232.97287197 3213.68213304

wfc3_tmin = [1943, 0, 1233]
wfc3_tmax = [3884, 2860, 3214]

#Spitzer Ch 2
#model, Tmin, Tmax =  zhang 1719.0 4257.78308959
#model, Tmin, Tmax =  hotspot_t 1380.0 3305.0e
#model, Tmin, Tmax =  spherical 930.545006798 3775.38559767
ch2_tmin = [1719., 1380, 930]
ch2_tmax = [4257, 3305, 3775]

#Spitzer Ch 1
#model, Tmin, Tmax =  zhang 1952.0 3672.60664016
#model, Tmin, Tmax =  hotspot_t 1427.0 3039.0
#model, Tmin, Tmax =  spherical 1277.04246551 3441.94916745
ch1_tmin = [1952., 1427., 1277.]
ch1_tmax = [3672., 3039., 3442.]


print "\\begin{deluxetable}{ccccccc}"
#\tabletypesize{\footnotesize} 
print "\\tablecolumns{7}"
print "\\tablewidth{0pt}:"
print "\\tablecaption{Model Comparison \label{table:models}}"
print "\\tablehead{"
print "\colhead{X} & \colhead{X} & \colhead{$T_\mathrm{min}$} & \colhead{$T_\mathrm{max}$} & \colhead{$\chi^2$} & \colhead{$\Delta$ AIC} & \colhead{$\Delta$ BIC}}"
print "\startdata"
print "\"Ch 1\" & \"m1\" & 0 & 0  & 1 & 1 & 1 \\\\"
print "\enddata"
print "\\vspace{-0.8cm}"
print "\\tablecomments{comments}"
print "\end{deluxetable}"

