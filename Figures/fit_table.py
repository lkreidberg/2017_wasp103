import numpy as np
from lc_fit import Model, LightCurveData
import pickle 

#get Spitzer AIC/BIC values from results.txt 

#zhang, hotspot, spherical, sincos
ch1_chi2 = np.array([0.9, 0.9, 0.9, 0.9])
ch1_aic = np.array([7131.0, 7107.4, 7096.7, 7047.7])
ch1_bic = np.array([7075.2, 7058.5, 7040.9, 7117.5])

ch2_chi2 = np.array([1.1, 1.1, 1.1, 1.1])
ch2_aic = np.array([8666.2, 8657.9, 8653.3, 8645.4])
ch2_bic = np.array([8722.0, 8706.8, 8709.1, 8729.1])


# bset fits for WFC3
wfc3_chi2 = []
wfc3_aic = []
wfc3_bic = []

models = ["zhang", "hotspot_t", "spherical", "sincos"]
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
#from get_mean_daynight_temp.py 
#zhang daymu, nightmu, min, max 2769.87230395 1976.73242901 1976.73242901 3953.15696275
#hotspot_t daymu, nightmu, min, max 2850.06092942 0.0 0.0 2878.84942366
#spherical daymu, nightmu, min, max 2635.76098954 1822.86068583 1227.05643006 3237.85099673
wfc3_tmin = np.array([1977, 0, 1227, 0])
wfc3_tmax = np.array([3953, 2879, 3237, 0])
wfc3_daymu = np.array([2769, 2879, 2636, 0])
wfc3_nightmu = np.array([1977, 0, 1822, 0])

#Spitzer Ch 2
#model, Tmin, Tmax =  zhang 1703.0 4276.80767467
#model, Tmin, Tmax =  hotspot_t 1361.0 3311.0
#model, Tmin, Tmax =  spherical 902.401726015 3786.1376002
ch2_tmin = np.array([1703, 1361, 902, 0])
ch2_tmax = np.array([4276, 3311, 3786, 0])
#from spitzer_maptemp.py
#zhang daymu, nightmu, min, max 2740.48761346 1711.12738853 1703.0 4287.03108111
#hotspot_t daymu, nightmu, min, max 3291.5 1361.0 1361.0 3311.0
#spherical daymu, nightmu, min, max 2919.98921277 1752.56989088 895.304735881 3786.3843194
ch2_daymu = np.array([2740, 3292, 2920, 0])
ch2_nightmu = np.array([1711, 1361, 1753, 0])

#Spitzer Ch 1
#model, Tmin, Tmax =  zhang 1955.0 3675.62794318
#model, Tmin, Tmax =  hotspot_t 1429.0 3041.0
#model, Tmin, Tmax =  spherical 1280.4447703 3444.26296074
ch1_tmin = np.array([1955, 1429, 1280, 0])
ch1_tmax = np.array([3675, 3041, 3444, 0])
#from spitzer_maptemp.py
#zhang daymu, nightmu, min, max 2656.25693185 1998.88866205 1955.0 3699.89072043
#hotspot_t daymu, nightmu, min, max 3024.88 1429.0 1429.0 3041.0
#spherical daymu, nightmu, min, max 2785.72145075 1935.79737607 1303.5292162 3424.67800668
ch1_daymu = np.array([2656, 3025, 2786, 0])
ch1_nightmu = np.array([1999, 1429, 1936, 0])


print "\\begin{deluxetable}{llllllll}"
#\tabletypesize{\footnotesize} 
print "\\tablecolumns{8}"
print "\\tablewidth{0pt}"
print "\\tablecaption{Model Comparison \label{table:models}}"
print "\\tablehead{"
print "\colhead{Data} & \colhead{Model} & \colhead{$T_\mathrm{min}$} & \colhead{$T_\mathrm{max}$} & \colhead{$\overline{T}_\mathrm{night}$} & \colhead{$\overline{T}_\mathrm{day}$} &\colhead{$\Delta_\mathrm{AIC}$} & \colhead{$\Delta_\mathrm{BIC}$}}"
print "\startdata"

# WFC 3
i = 2
print "WFC3 & Sph. Harmonics &", str(wfc3_tmin[i]), "&", str(wfc3_tmax[i]), "&",  str(wfc3_nightmu[i]), "&", str(wfc3_daymu[i]), "&",np.round(wfc3_aic[i] - wfc3_aic.min(), decimals = 1), "&", np.round(wfc3_bic[i] - wfc3_bic.min(), decimals = 1), "\\\\"
i = 0
print "\, & Kinematic &", str(wfc3_tmin[i]), "&", str(wfc3_tmax[i]), "&", str(wfc3_nightmu[i]), "&", str(wfc3_daymu[i]), "&",np.round(wfc3_aic[i] - wfc3_aic.min(), decimals = 1), "&", np.round(wfc3_bic[i] - wfc3_bic.min(), decimals = 1), "\\\\"
i = 1
print "\, & Two Temp. &", str(wfc3_tmin[i]), "&", str(wfc3_tmax[i]), "&",  str(wfc3_nightmu[i]), "&", str(wfc3_daymu[i]), "&",np.round(wfc3_aic[i] - wfc3_aic.min(), decimals = 1), "&", np.round(wfc3_bic[i] - wfc3_bic.min(), decimals = 1), "\\\\"
i = 3
print "\, & Sinusoid &", "--", "&", "--", "&",  "--", "&", "--", "&",np.round(wfc3_aic[i] - wfc3_aic.min(), decimals = 1), "&", np.round(wfc3_bic[i] - wfc3_bic.min(), decimals = 1), "\\\\"

# Spitzer Channel 1
i = 2
print "Ch 1 & Sph. Harmonics &", str(ch1_tmin[i]), "&", str(ch1_tmax[i]), "&", str(ch1_nightmu[i]), "&", str(ch1_daymu[i]), "&",str(ch1_aic[i] - ch1_aic.min()), "&", str(ch1_bic[i] - ch1_bic.min()), "\\\\"
i = 0
print "\, & Kinematic &", str(ch1_tmin[i]), "&", str(ch1_tmax[i]), "&", str(ch1_nightmu[i]), "&", str(ch1_daymu[i]), "&",str(ch1_aic[i] - ch1_aic.min()), "&", str(ch1_bic[i] - ch1_bic.min()), "\\\\"
i = 1
print "\, & Two Temp. &", str(ch1_tmin[i]), "&", str(ch1_tmax[i]), "&", str(ch1_nightmu[i]), "&", str(ch1_daymu[i]), "&",str(ch1_aic[i] - ch1_aic.min()), "&", str(ch1_bic[i] - ch1_bic.min()), "\\\\"
i = 3
print "\, & Sinusoid &", "--", "&", "--", "&", "--", "&", "--", "&", np.round(ch1_aic[i] - ch1_aic.min(), decimals = 1), "&", np.round(ch1_bic[i] - ch1_bic.min(), decimals = 1), "\\\\"


#Spitzer Channel 2
i = 2
print "Ch 2 & Sph. Harmonics &", str(ch2_tmin[i]), "&", str(ch2_tmax[i]), "&", str(ch2_nightmu[i]), "&", str(ch2_daymu[i]), "&",str(ch2_aic[i] - ch2_aic.min()), "&", str(ch2_bic[i] - ch2_bic.min()), "\\\\"
i = 0
print "\, & Kinematic &", str(ch2_tmin[i]), "&", str(ch2_tmax[i]), "&",  str(ch2_nightmu[i]), "&", str(ch2_daymu[i]), "&",str(ch2_aic[i] - ch2_aic.min()), "&", str(ch2_bic[i] - ch2_bic.min()), "\\\\"
i = 1
print "\, & Two Temp. &", str(ch2_tmin[i]), "&", str(ch2_tmax[i]), "&", str(ch2_nightmu[i]), "&", str(ch2_daymu[i]), "&", str(ch2_aic[i] - ch2_aic.min()), "&", str(ch2_bic[i] - ch2_bic.min()), "\\\\"
i = 3
print "\, & Sinusoid &", "--", "&", "--", "&", "--", "&", "--", "&", np.round(ch2_aic[i] - ch2_aic.min(), decimals = 1), "&", np.round(ch2_bic[i] - ch2_bic.min(), decimals = 1), "\\\\"

print "\enddata"
print "\\vspace{-0.8cm}"
print "\\tablecomments{comments}"
print "\end{deluxetable}"

