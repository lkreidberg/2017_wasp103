import numpy as np
import matplotlib.pyplot as plt


def blackbody(l,T): 
	h = 6.62607e-34   #J/s
        c = 2.997925e8       #m/s
        k = 1.38065e-23     #J/K
        return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T))-1))                 #[W/sr/m^3]

w = np.linspace(1.1, 1.7, 100)
w = np.linspace(5, 10, 100)
star = np.genfromtxt("wasp103_sed_fluxes.out")
star_bb = np.interp(w, star[:,0], star[:,1])*1.e24/(w*np.pi*4.)
plt.plot(w, star_bb)

Tstar = 6110.
star_bb = blackbody(w*1.0e-6, Tstar)
plt.plot(w, star_bb)

#plt.gca().set_yscale('log')

plt.show()
