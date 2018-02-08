import numpy as np
from lc_fit import Model, LightCurveData
import pickle 

#get Spitzer AIC/BIC values from results.txt 

#zhang, hotspot, spherical, sincos
ch1_chi2 = np.array([0.9, 0.9, 0.9, 0.9])
ch1_aic = np.array([7075.1, 7058.5, 7040.9, 7044.7])
ch1_bic = np.array([7131.1, 7107.4, 7096.7, 7121.4])


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
#from spitzer_maptemp.py
#zhang daymu, nightmu, min, max 2544.42433961 1621.25879891 1614.0 3931.36778946
#hotspot_t daymu, nightmu, min, max 3222.03 1344.0 1344.0 3241.0
#spherical daymu, nightmu, min, max 2864.26860292 1728.83922558 888.59900134 3713.55482549
ch2_daymu = np.array([2544, 3241, 2864, 0])
ch2_nightmu = np.array([1621, 1344, 1729, 0])
ch2_tmin = np.array([1614, 1344, 888, 0])
ch2_tmax = np.array([3931, 3241, 3714, 0])

#Spitzer Ch 1
#from spitzer_maptemp.py
#models = ["zhang", "hotspot_t", "spherical", "sincos"]
#zhang daymu, nightmu, min, max 2614.33231628 1974.7502736 1932.0 3629.80700341
#hotspot_t daymu, nightmu, min, max 2974.28 1418.0 1418.0 2990.0
#spherical daymu, nightmu, min, max 2740.6544324 1912.86267538 1268.93864177 3391.17897337
ch1_daymu = np.array([2614, 2990, 2741, 0])
ch1_nightmu = np.array([1975, 1418, 1912, 0])
ch1_tmin = np.array([1932, 1418, 1269, 0])
ch1_tmax = np.array([3630, 2990, 3391, 0])


print "\\begin{deluxetable}{lllllll}"
#\tabletypesize{\footnotesize} 
print "\\tablecolumns{7}"
print "\\tablewidth{0pt}"
print "\\tablecaption{Model Comparison \label{table:models}}"
print "\\tablehead{"
print "\colhead{Data} & \colhead{Model} & \colhead{$T_\mathrm{min}$} & \colhead{$T_\mathrm{max}$} & \colhead{$\overline{T}_\mathrm{night}$} & \colhead{$\overline{T}_\mathrm{day}$} & \colhead{$\Delta_\mathrm{BIC}$}}"
print "\startdata"

# WFC 3
i = 2
print "WFC3 & Sph. Harmonics &", str(wfc3_tmin[i]), "&", str(wfc3_tmax[i]), "&",  str(wfc3_nightmu[i]), "&", str(wfc3_daymu[i]), "&", np.round(wfc3_bic[i] - wfc3_bic.min(), decimals = 1), "\\\\"
i = 0
print "\, & Kinematic &", str(wfc3_tmin[i]), "&", str(wfc3_tmax[i]), "&", str(wfc3_nightmu[i]), "&", str(wfc3_daymu[i]), "&", np.round(wfc3_bic[i] - wfc3_bic.min(), decimals = 1), "\\\\"
i = 1
print "\, & Two Temp. &", str(wfc3_tmin[i]), "&", str(wfc3_tmax[i]), "&",  str(wfc3_nightmu[i]), "&", str(wfc3_daymu[i]), "&", np.round(wfc3_bic[i] - wfc3_bic.min(), decimals = 1), "\\\\"
i = 3
print "\, & Sinusoid &", "--", "&", "--", "&",  "--", "&", "--", "&", np.round(wfc3_bic[i] - wfc3_bic.min(), decimals = 1), "\\\\"

# Spitzer Channel 1
i = 2
print "Ch 1 & Sph. Harmonics &", str(ch1_tmin[i]), "&", str(ch1_tmax[i]), "&", str(ch1_nightmu[i]), "&", str(ch1_daymu[i]), "&", str(ch1_bic[i] - ch1_bic.min()), "\\\\"
i = 0
print "\, & Kinematic &", str(ch1_tmin[i]), "&", str(ch1_tmax[i]), "&", str(ch1_nightmu[i]), "&", str(ch1_daymu[i]), "&", str(ch1_bic[i] - ch1_bic.min()), "\\\\"
i = 1
print "\, & Two Temp. &", str(ch1_tmin[i]), "&", str(ch1_tmax[i]), "&", str(ch1_nightmu[i]), "&", str(ch1_daymu[i]), "&", str(ch1_bic[i] - ch1_bic.min()), "\\\\"
i = 3
print "\, & Sinusoid &", "--", "&", "--", "&", "--", "&", "--", "&", np.round(ch1_bic[i] - ch1_bic.min(), decimals = 1), "\\\\"


#Spitzer Channel 2
i = 2
print "Ch 2 & Sph. Harmonics &", str(ch2_tmin[i]), "&", str(ch2_tmax[i]), "&", str(ch2_nightmu[i]), "&", str(ch2_daymu[i]), "&", str(ch2_bic[i] - ch2_bic.min()), "\\\\"
i = 0
print "\, & Kinematic &", str(ch2_tmin[i]), "&", str(ch2_tmax[i]), "&",  str(ch2_nightmu[i]), "&", str(ch2_daymu[i]), "&", str(ch2_bic[i] - ch2_bic.min()), "\\\\"
i = 1
print "\, & Two Temp. &", str(ch2_tmin[i]), "&", str(ch2_tmax[i]), "&", str(ch2_nightmu[i]), "&", str(ch2_daymu[i]), "&",  str(ch2_bic[i] - ch2_bic.min()), "\\\\"
i = 3
print "\, & Sinusoid &", "--", "&", "--", "&", "--", "&", "--", "&",  np.round(ch2_bic[i] - ch2_bic.min(), decimals = 1), "\\\\"

print "\enddata"
print "\\vspace{-0.8cm}"
print "\\tablecomments{comments}"
print "\end{deluxetable}"

