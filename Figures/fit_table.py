import numpy as np
from lc_fit import Model, LightCurveData
import pickle 
#zhang, hotspot, spherical
ch1_chi2 = np.array([0.9, 0.9, 0.9])
ch1_aic = np.array([7131.0, 7107.4, 7096.7])
ch1_bic = np.array([7075.2, 7058.5, 7040.9])

ch2_chi2 = np.array([1.1, 1.1, 1.1])
ch2_aic = np.array([8345.95, 8355.5, 8331.6])
ch2_bic = np.array([8401.61, 8384.2, 8387.3])

# bset fits for WFC3
wfc3_chi2 = []
wfc3_aic = []
wfc3_bic = []

models = ["zhang", "hotspot_t", "spherical"]
path = "WFC3_best_fits/"
for i, model in enumerate(models):
	p = pickle.load(open(path+"bestfit_"+model+".pic", "rb")) 

	d, m, par = p[0], p[1], p[2]            #stores data,  model, and best fit parameters into d, m, & par
	wfc3_chi2.append(m.chi2red)
	wfc3_bic.append(m.chi2 + d.nfree_param*np.log(d.npoints))
	wfc3_aic.append(m.chi2 + 2.*d.nfree_param)
#	print m.chi2red, bic, aic

wfc3_chi2 = np.array(wfc3_chi2)
wfc3_aic = np.array(wfc3_aic)
wfc3_bic = np.array(wfc3_bic)

#Temperatures:
#HST:
#model, Tmin, Tmax =  zhang 1979.6612603 3958.07678595
#model, Tmin, Tmax =  hotspot_t 0.0 2887.11857794
#model, Tmin, Tmax =  spherical 1213.50410565 3251.37484203
wfc3_tmin = np.array([1980, 0, 1214])
wfc3_tmax = np.array([3958, 2887, 3251])

#Spitzer Ch 2
#model, Tmin, Tmax =  zhang 1703.0 4276.80767467
#model, Tmin, Tmax =  hotspot_t 1361.0 3311.0
#model, Tmin, Tmax =  spherical 902.401726015 3786.1376002
ch2_tmin = np.array([1703, 1361, 902])
ch2_tmax = np.array([4276, 3311, 3786])

#Spitzer Ch 1
#model, Tmin, Tmax =  zhang 1955.0 3675.62794318
#model, Tmin, Tmax =  hotspot_t 1429.0 3041.0
#model, Tmin, Tmax =  spherical 1280.4447703 3444.26296074
ch1_tmin = np.array([1955, 1429, 1280])
ch1_tmax = np.array([3675, 3041, 3444])


print "\\begin{deluxetable}{lllllll}"
#\tabletypesize{\footnotesize} 
print "\\tablecolumns{7}"
print "\\tablewidth{0pt}:"
print "\\tablecaption{Model Comparison \label{table:models}}"
print "\\tablehead{"
print "\colhead{Data} & \colhead{Model} & \colhead{$T_\mathrm{min}$} & \colhead{$T_\mathrm{max}$} & \colhead{$\chi^2$} & \colhead{$\Delta_\mathrm{AIC}$} & \colhead{$\Delta_\mathrm{BIC}$}}"
print "\startdata"

# WFC 3
i = 2
print "WFC3 & Sph. Harmonics &", str(wfc3_tmin[i]), "&", str(wfc3_tmax[i]), "&", np.round(wfc3_chi2[i], decimals = 1), "&", np.round(wfc3_aic[i] - wfc3_aic.min(), decimals = 1), "&", np.round(wfc3_bic[i] - wfc3_bic.min(), decimals = 1), "\\\\"
i = 0
print "\, & Kinematic &", str(wfc3_tmin[i]), "&", str(wfc3_tmax[i]), "&", np.round(wfc3_chi2[i], decimals = 1), "&", np.round(wfc3_aic[i] - wfc3_aic.min(), decimals = 1), "&", np.round(wfc3_bic[i] - wfc3_bic.min(), decimals = 1), "\\\\"
i = 1
print "\, & Two Temp. &", str(wfc3_tmin[i]), "&", str(wfc3_tmax[i]), "&", np.round(wfc3_chi2[i], decimals = 1), "&", np.round(wfc3_aic[i] - wfc3_aic.min(), decimals = 1), "&", np.round(wfc3_bic[i] - wfc3_bic.min(), decimals = 1), "\\\\"

# Spitzer Channel 1
i = 2
print "Ch 1 & Sph. Harmonics &", str(ch1_tmin[i]), "&", str(ch1_tmax[i]), "&", str(ch1_chi2[i]), "&", str(ch1_aic[i] - ch1_aic.min()), "&", str(ch1_bic[i] - ch1_bic.min()), "\\\\"
i = 0
print "\, & Kinematic &", str(ch1_tmin[i]), "&", str(ch1_tmax[i]), "&", str(ch1_chi2[i]), "&", str(ch1_aic[i] - ch1_aic.min()), "&", str(ch1_bic[i] - ch1_bic.min()), "\\\\"
i = 1
print "\, & Two Temp. &", str(ch1_tmin[i]), "&", str(ch1_tmax[i]), "&", str(ch1_chi2[i]), "&", str(ch1_aic[i] - ch1_aic.min()), "&", str(ch1_bic[i] - ch1_bic.min()), "\\\\"


#Spitzer Channel 2
i = 2
print "Ch 2 & Sph. Harmonics &", str(ch2_tmin[i]), "&", str(ch2_tmax[i]), "&", str(ch2_chi2[i]), "&", str(ch2_aic[i] - ch2_aic.min()), "&", str(ch2_bic[i] - ch2_bic.min()), "\\\\"
i = 0
print "\, & Kinematic &", str(ch2_tmin[i]), "&", str(ch2_tmax[i]), "&", str(ch2_chi2[i]), "&", str(ch2_aic[i] - ch2_aic.min()), "&", str(ch2_bic[i] - ch2_bic.min()), "\\\\"
i = 1
print "\, & Two Temp. &", str(ch2_tmin[i]), "&", str(ch2_tmax[i]), "&", str(ch2_chi2[i]), "&", str(ch2_aic[i] - ch2_aic.min()), "&", str(ch2_bic[i] - ch2_bic.min()), "\\\\"

print "\enddata"
print "\\vspace{-0.8cm}"
print "\\tablecomments{comments}"
print "\end{deluxetable}"

