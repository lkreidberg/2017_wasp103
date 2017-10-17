import numpy as np
import matplotlib.pyplot as plt

def blackbody(l,T): 
	h = 6.626e-27 	#ergs*s
	c = 3.e10	#cm/s
	kb = 1.38e-16	#erg/K
	return 2*h*c**2/(l**5*(np.exp(h*c/(l*kb*T))-1))

d = np.genfromtxt("wasp103_sed_fluxes.out")
w, f  = d[:,0], d[:,1]*np.pi*4.		#multiplies flux by steradians

bb = blackbody(w/1.e4, 6110.)

print f[903], bb[903]

#plt.plot(d[:,0], d[:,1])

#plt.xlim(1, 5)


