import numpy as np
import matplotlib.pyplot as plt

def blackbody(l,T): 
	h = 6.62607e-34   #J/s
        c = 2.997925e8       #m/s
        k = 1.38065e-23     #J/K
        return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T))-1))                 #[W/sr/m^3]

d =  np.genfromtxt("W103-6110K-1.22Ms-1.22Rs.dat")

outfile = open("w103_spectrum.txt", "w")
for i in range(len(d)): print>>outfile, d[i,0]*1e-6, d[i,1]*22423.*1.e6/(np.pi)
outfile.close()


d = np.genfromtxt("w103_spectrum.txt")


Ts = 6110.
plt.plot(d[:,0], d[:,1])
plt.plot(d[:,0], blackbody(d[:,0], Ts))
plt.xlim(0,3e-6)

plt.show()
