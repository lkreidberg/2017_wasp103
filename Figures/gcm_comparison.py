import matplotlib.pyplot as plt
import numpy as np
import os, glob
from lc_fit import Model, LightCurveData
import pickle

gcm_path = "./Vivien_models2"
gcms = glob.glob(os.path.join(gcm_path, "PC*PT.dat"))

labels = ["NoTiO", "Drag1", "Drag2", "Drag3", "Met0.5", "TiO"]


def proj_area(phi, inc):
        R0 = 1.07e8             #planet radius in m
        p = 0.00117             # planet-to-star mass ratio
        r = 2.9695e9            #orbital separation in m
        kn = 0.653              #from table B5 of Leconte et al. 2011
        n = 1.005
        qn = kn*(1.- n/5.)      #n 


        alpha1 = 2.5*qn*R0**3/(p*r**3)
        alpha2 = -1.25*qn*R0**3/(p*r**3)
        alpha3 = alpha2

        a1 = R0*(1+alpha1)
        a2 = R0*(1+alpha2)
        a3 = R0*(1+alpha3)

        return  np.sqrt(a3**2*np.sin(inc)**2*(a1**2*np.sin(phi)**2+a2**2*np.cos(phi)**2)+ a1**2*a2**2*np.cos(inc)**2)/(a2*a3)


inc = 1.523

plt.subplot(211)
#plot gcms
for i, g in enumerate(gcms):
    print g,
    model = np.genfromtxt(g, delimiter = ',') 
    phi = model[:,0]/360.*2.*np.pi
    area = proj_area(phi, inc)
    plt.plot(model[:,0]/360. +0.5, model[:,5]*1e3, label = labels[i])
    print (np.max(model[:,5]) - np.min(model[:,1]))*1e3
    


f = "./WFC3_best_fits/old_white_fits/bestfit_spherical.pic"
p = pickle.load(open(f, 'rb'))
d, m, par = p[0], p[1], p[2]		#stores data,  model, and best fit parameters into d, m, & par
dilution = d.dilution + 1.

ind = d.err < 9.0e7			#indices of outliers

err = d.err[ind]/d.flux[ind]		#normalized per point uncertainty 
phase = m.phase[ind]
data_corr = m.data_corr[ind]
bestfit = m.bestfit[ind]


ind = np.argsort(phase)
err, phase, data_corr, bestfit = err[ind], phase[ind], data_corr[ind], bestfit[ind] #sorts by phase

plt.plot(phase, (bestfit-1.)*1e3*dilution, color = 'k', label = 'best fit')

plt.legend(loc = 'upper right')
plt.ylim(0., 2)
plt.ylabel("fp/fs")
plt.xlabel("Phase")
plt.title("No ellipsoidal correction")


plt.subplot(212)
#plot gcms
for i, g in enumerate(gcms):
    print g,
    model = np.genfromtxt(g, delimiter = ',') 
    phi = model[:,0]/360.*2.*np.pi
    area = proj_area(phi, inc)
    plt.plot(model[:,0]/360. +0.5, model[:,5]*1e3*area, label = labels[i])
    print (np.max(model[:,5]) - np.min(model[:,1]))*1e3


f = "./WFC3_best_fits/old_white_fits/bestfit_spherical.pic"
p = pickle.load(open(f, 'rb'))
d, m, par = p[0], p[1], p[2]		#stores data,  model, and best fit parameters into d, m, & par
dilution = d.dilution + 1.

ind = d.err < 9.0e7			#indices of outliers

err = d.err[ind]/d.flux[ind]		#normalized per point uncertainty 
phase = m.phase[ind]
data_corr = m.data_corr[ind]
bestfit = m.bestfit[ind]


ind = np.argsort(phase)
err, phase, data_corr, bestfit = err[ind], phase[ind], data_corr[ind], bestfit[ind] #sorts by phase

plt.plot(phase, (bestfit-1.)*1e3*dilution, color = 'k', label = 'best fit')
plt.ylim(0., 2)
plt.ylabel("fp/fs")
plt.xlabel("Phase")
plt.legend(loc = 'upper right')
plt.title("Ellipsoidal correction")

plt.tight_layout()
plt.savefig("gcm_comparison.pdf")
plt.show()
