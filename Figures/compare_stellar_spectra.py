import numpy as np
import matplotlib.pyplot as plt

#Vivien spectrum - col = wavelength (microns), flux at 1 AU, in W/m2/micron, flux at stellar surface, erg/cm2/sec/Hz 
d = np.genfromtxt("W103-6110K-1.22Ms-1.22Rs.dat")
d[:,1] = d[:,1]*d[:,0]
plt.plot(d[:,0], d[:,1], label = "Vivien")
print d[:,1].max()

#Keivan spectrum, micron, [ergs/s/cm^2]  
d = np.genfromtxt("wasp103_sed_fluxes.out")
plt.plot(d[:,0], d[:,1], label = "Keivan")
plt.xlim(0,5)
plt.legend(loc= 'upper right')
print d[:,1].max()

plt.show()
