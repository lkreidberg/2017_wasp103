import matplotlib.pyplot as plt
import numpy as np
import os, glob
from lc_fit import Model, LightCurveData
import pickle
import seaborn as sns
from scipy.interpolate import interp1d

sns.set_context("paper")
sns.set_style("white")
sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})

gcm_path = "./Vivien_models2"
gcms = glob.glob(os.path.join(gcm_path, "PC*PT.dat"))

labels = ["$\\tau_\mathrm{drag} = 10^4$ s", "$\\tau_\mathrm{drag} = 10^3$ s", "[Fe/H] = 0.5", "nominal GCM"]

plt.figure(figsize = (4, 2.5))

#colors = ["#be0119", "#f97306", "#04d9ff", "#0165fc"]
colors = ["#be0119", "#f97306", "#6dedfd", "#0165fc"]
ls = [":", "--", "-.", "-"]


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

#plot gcms
for i, g in enumerate(gcms):
    print g,
    model = np.genfromtxt(g, delimiter = ',') 
    phi = model[:,0]/360.*2.*np.pi
    area = proj_area(phi, inc)
    plt.plot(model[:,0]/360. +0.5, model[:,5]*1e3*area, label = labels[i], color = colors[i], linestyle = ls[i])
    plt.plot(model[:,0]/360. +1.5, model[:,5]*1e3*area, color = colors[i], linestyle = ls[i])

    #print "model quadrature", np.interp(0.25, model[:,0]/360. +0.5, model[:,5]*1e3*area) 

    #for col in [5,6,7]:
    #col = 5
    #col =6
    col = 7
    x, y = model[:,0]/360. + 0.5, model[:,col]
    f = interp1d(x, y, kind = 'cubic')
    xnew = np.linspace(0, 1, 1000)
    ynew = f(xnew)
    ind = ynew == np.max(ynew)
    offset = xnew[ind]
    print ((offset - 0.5)*(-360.))[0], (np.max(ynew) - np.min(ynew))/np.max(ynew),
    print ""



#f = "./WFC3_best_fits/old_white_fits/bestfit_spherical.pic"
f = "./WFC3_best_fits//bestfit_sincos.pic"
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

plt.plot(phase, (bestfit-1.)*1e3*dilution, color = '0.2', label = 'best fit')
plt.plot(phase + 1., (bestfit-1.)*1e3*dilution, color = '0.2') 

print "quadrature best fit", np.interp(0.25, phase, (bestfit - 1.)*1e3*dilution)

plt.xlim(0, 1.1)
plt.ylim(-0.1, 1.8)
plt.ylabel("Planet-to-star flux (ppt)")
plt.xlabel("Phase")

ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='upper right')

plt.tight_layout()
plt.savefig("gcm_comparison.pdf")
#plt.show()
