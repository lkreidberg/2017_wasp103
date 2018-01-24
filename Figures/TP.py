import numpy as np
import matplotlib.pyplot as plt
import pickle
#import seaborn as sns
from pylab import *

#sns.set_context("talk")
#sns.set_style("white")
#sns.set_style("ticks", {"xtick.direction":"in", "ytick.direction":"in"})

rc('font',**{'family':'sans-serif','sans-serif':['Arial']})


plt.figure(figsize=(4,6))

color = "#d9544d"

#Dayside
plt.subplot(211)

data_wl, data, data_err,  best_fit_binned,  wl_hi,  y_hi_best, spec_arr, Tarr, P, samples = pickle.load(open("Mike_models/WASP-103b_grid_DAYSIDE_output.pic", "rb"))

levels = [7, 14, 20, 27]
for l in levels: print P[l], np.percentile(Tarr[:,l], np.array([16., 50., 84.]))

for i in range(99):
    plt.plot(Tarr[i,:], P, color = color)

plt.plot(Tarr[99,:], P, color = color, label = '1-D models')

d = np.genfromtxt("GCM_From_Vivien/PTprofiles-WASP-103b-TiO-fix-3-Drag3-NEW-OPA-nit-1036800.dat", delimiter = ",")
#Layer number, Latitude, Longitude, Pressure(bars), Temperature(K), -1,-1,-1,-1

lat, lon = d[:,1], d[:,2]
latmed = np.percentile(lat, 50., interpolation='nearest')
lonmed = np.percentile(lon, 50., interpolation='nearest')
ind = (lat==latmed)&(lon == lonmed)

plt.plot(d[ind,4], d[ind,3], color = 'k', label = 'GCM')
plt.legend(loc = 'upper left')

plt.gca().set_yscale('log')
plt.ylim(1e1, 2e-4)
plt.xlim(2200, 3800)

plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (bar)')
plt.title("Dayside")


#Nightside
plt.subplot(212)
data_wl, data, data_err,  best_fit_binned,  wl_hi,  y_hi_best, spec_arr, Tarr, P, samples = pickle.load(open("Mike_models/WASP-103b_grid_DAYSIDE_output.pic", "rb"))

levels = [7, 14, 20, 27]
for l in levels: print P[l], np.percentile(Tarr[:,l], np.array([16., 50., 84.]))

"""for i in range(99):
    plt.plot(Tarr[i,:], P, color = color)

plt.plot(Tarr[99,:], P, color = color, label = '1-D models')"""

d = np.genfromtxt("GCM_From_Vivien/PTprofiles-WASP-103b-TiO-fix-3-Drag3-NEW-OPA-nit-1036800.dat", delimiter = ",")
#Layer number, Latitude, Longitude, Pressure(bars), Temperature(K), -1,-1,-1,-1

lat, lon = d[:,1], d[:,2]
latmed = np.percentile(lat, 50., interpolation='nearest')
lonmed = np.percentile(lon, 50., interpolation='nearest')
ind = (lat==latmed)&(lon == 81.562)

plt.plot(d[ind,4], d[ind,3], color = 'k', label = 'GCM')
plt.legend(loc = 'upper left')

plt.gca().set_yscale('log')
plt.ylim(1e1, 2e-4)
#plt.xlim(2200, 3800)

plt.title("Nightside")
plt.xlabel('Temperature (K)')
plt.ylabel('Pressure (bar)')





plt.tight_layout()
plt.savefig("TP.pdf")
#plt.show()
