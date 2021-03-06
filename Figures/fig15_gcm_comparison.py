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
colors = ["#be0119", "#f97306", "#6dedfd", "#0165fc"]
ls = [":", "--", "-.", "-"]

plt.figure(figsize = (4, 6))



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


plt.subplot(311)
inc = 1.523
correction_factor = 1.099  # = (0.1146/0.1093)**2, ratio of correct transit depth to what Vivien used (which was from the discovery paper and not corrected for stellar contamination)

#plot gcms
for i, g in enumerate(gcms):
    print g,
    model = np.genfromtxt(g, delimiter = ',') 
    phi = model[:,0]/360.*2.*np.pi
    area = proj_area(phi, inc)
    plt.plot(model[:,0]/360. +0.5, model[:,5]*1e3*area*correction_factor, color = colors[i], linestyle = ls[i])
    plt.plot(model[:,0]/360. +1.5, model[:,5]*1e3*area*correction_factor, color = colors[i], linestyle = ls[i])

    #get amlitude and offset for different gcms
    col = 7							#5,6,7 = WFC3, Ch1, Ch2
    x, y = model[:,0]/360. + 0.5, model[:,col]
    f = interp1d(x, y, kind = 'cubic')
    xnew = np.linspace(0, 1, 1000)
    ynew = f(xnew)
    ind = ynew == np.max(ynew)
    offset = xnew[ind]
    print ((offset - 0.5)*(-360.))[0], (np.max(ynew) - np.min(ynew))/np.max(ynew),
    print ""

labels = ["$\\tau_\mathrm{drag} = 10^4$ s", "$\\tau_\mathrm{drag} = 10^3$ s", "[Fe/H] = 0.5", "nominal GCM"]
colors = ["#be0119", "#f97306", "#6dedfd", "#0165fc"]

x, y = np.linspace(100, 101, 10), np.linspace(100, 101, 10)
l1, =	plt.plot(x,y, colors[3], linestyle = ls[3], label = labels[3])
l2, =	plt.plot(x,y, colors[2], linestyle = ls[2], label = labels[2])
l3, =	plt.plot(x,y, colors[1], linestyle = ls[1], label = labels[1])
l4, =	plt.plot(x,y, colors[0], linestyle = ls[0], label = labels[0])


f = "./WFC3_best_fits/bestfit_spherical.pic"
#f = "./WFC3_best_fits//bestfit_sincos.pic"
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

#l5, = plt.plot(phase, (bestfit-1.)*1e3*dilution, color = '0.2', label = 'best fit sine curve')
l5, = plt.plot(phase, (bestfit-1.)*1e3*dilution, color = '0.2', label = 'best fit model')
plt.plot(phase + 1., (bestfit-1.)*1e3*dilution, color = '0.2') 

phasebins = np.linspace(0., 1., 30)
nbins = len(phasebins) - 1

bindata = np.zeros(nbins)
binsigma = np.zeros(nbins)
binbestfit = np.zeros(nbins)
bin_average = np.zeros(nbins)

for j in range(1, len(phasebins)):
        ind = (phase >= phasebins[j-1]) & (phase < phasebins[j])
        if sum(ind)>0:
                bindata[j-1] = sum(data_corr[ind]/err[ind]**2)/sum(1/err[ind]**2) 
                binsigma[j-1] = np.sqrt(1/sum(1/err[ind]**2))
                binbestfit[j-1]=  np.mean(bestfit[ind])
                bin_average[j-1] = (phasebins[j-1]+phasebins[j])/2.

plt.errorbar(bin_average, (bindata-1.)*1e3*dilution, yerr = binsigma*1e3, fmt = '.k')

plt.xlim(0, 1.1)
plt.ylim(-0.1, 2.3)
plt.ylabel("Planet-to-star flux ($\\times10^{-3}$)")
plt.xlabel("Phase")

first_legend = plt.legend(handles=[l1, l2, l3, l4], loc='upper right', frameon=True)
ax = plt.gca().add_artist(first_legend)

#second_legend = plt.legend(handles=[l5], loc=(0.26, 0.834), frameon=True)
second_legend = plt.legend(handles=[l5], loc='upper left', frameon=True)
ax = plt.gca().add_artist(second_legend)

plt.gca().text(0.05, 1.6, 'HST/WFC3', fontsize=10)

plt.subplot(312)

#Spitzer Ch 1
for i, g in enumerate(gcms):
    print g,
    model = np.genfromtxt(g, delimiter = ',') 
    phi = model[:,0]/360.*2.*np.pi
    area = proj_area(phi, inc)
    plt.plot(model[:,0]/360. +0.5, model[:,6]*1e3*area*correction_factor, label = labels[i], color = colors[i], linestyle = ls[i])
    plt.plot(model[:,0]/360. +1.5, model[:,6]*1e3*area*correction_factor, color = colors[i], linestyle = ls[i])

plt.xlim(0, 1.1)
plt.ylim(-0.6, 5.5)
plt.ylabel("Planet-to-star flux ($\\times10^{-3}$)")
plt.xlabel("Phase")

plt.gca().text(0.05, 4., 'Spitzer\nCh 1', fontsize=10)

plt.subplot(313)

#Spitzer Ch 2
for i, g in enumerate(gcms):
    print g,
    model = np.genfromtxt(g, delimiter = ',') 
    phi = model[:,0]/360.*2.*np.pi
    area = proj_area(phi, inc)
    plt.plot(model[:,0]/360. +0.5, model[:,7]*1e3*area*correction_factor, label = labels[i], color = colors[i], linestyle = ls[i])
    plt.plot(model[:,0]/360. +1.5, model[:,7]*1e3*area*correction_factor, color = colors[i], linestyle = ls[i])

plt.xlim(0, 1.1)
plt.ylim(-0.6, 6.3)
plt.ylabel("Planet-to-star flux ($\\times10^{-3}$)")
plt.xlabel("Phase")

plt.gca().text(0.05, 4.5, 'Spitzer\nCh 2', fontsize=10)


#plot spitzer best fits and data

dilution = np.array([0.1712, 0.1587])
l1 = np.array([3.15e-6, 4.0e-6])
l2 = np.array([3.95e-6, 5.0e-6])
#bestfits = ["Ch1_best_fits/2017-09-15_11:33-sincos2/bestfit.pic", "Ch2_best_fits/2017-09-15_12:02-sincos2/bestfit.pic"]
bestfits = ["Ch1_best_fits/2018-02-07_14:28-spherical/bestfit.pic", "Ch2_best_fits/2018-02-07_14:07-spherical/bestfit.pic"]

N = 2000

for ii, f in enumerate(bestfits):
        plt.subplot(312 + ii)
	p = pickle.load(open(f, 'rb'))
	data_corr = p[1] 
	err = p[2] 
	bestfit = p[3]
	phase = p[0] 

        data_corr, err, bestfit, phase = data_corr[N::], err[N::], bestfit[N::], phase[N::]


	ind = phase > 1.0
	phase[ind] -= 1.0

        ind = (phase>0.49)&(phase<0.51)
        offset = np.mean(bestfit[ind])
        print "offset", offset
        bestfit -= offset
        data_corr -= offset
        
	bindata = np.zeros(nbins)
	binsigma = np.zeros(nbins)
	binbestfit = np.zeros(nbins)
	bin_average = np.zeros(nbins)
	
	for j in range(1, len(phasebins)):
		ind = (phase >= phasebins[j-1]) & (phase < phasebins[j])
		if sum(ind)>0:
			bindata[j-1] = sum(data_corr[ind]/err[ind]**2)/sum(1/err[ind]**2) 
			binsigma[j-1] = np.sqrt(1/sum(1/err[ind]**2))
			binbestfit[j-1]=  np.mean(bestfit[ind])
			bin_average[j-1] = (phasebins[j-1]+phasebins[j])/2.

	plt.errorbar(bin_average, bindata*1e3*(1.+dilution[ii]), yerr = binsigma*1e3, marker = '.', linestyle = 'None', color = '0.2')
	plt.plot(phase, bestfit*1e3*(1.+dilution[ii]), color = '0.2')





plt.tight_layout()
plt.savefig("fig15.pdf")
#plt.show()
