import numpy as np
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel, convolve

def blackbody(l,T): 
	h = 6.62607e-34   #J/s
        c = 2.997925e8       #m/s
        k = 1.38065e-23     #J/K
        return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T))-1))                 #[W/sr/m^3]

g = Gaussian1DKernel(stddev=150)

#Vivien's stellar spectrum
d = np.genfromtxt("W103-6110K-1.22Ms-1.22Rs.dat")       #W/m2/micron (column 1)
plt.plot(d[:,0], convolve(d[:,1]*22423., g))             #multiply by (1 astronomical unit/1.436 solar radii)**2

#blackbody model
w = np.linspace(0.1, 10., 1000)*1.e-6
bb  = blackbody(w, 6110.)

plt.plot(w*1.e6, np.pi*bb/1.e6)             #divide by 1e6 to get in units of W/m2/*micron*, multiple by pi for steradians

plt.xlim(0, 5)
plt.show()
