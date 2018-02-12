import numpy as np
import spiderman_spherical
import spiderman
import matplotlib.pyplot as plt

def quantile(x, q): return np.percentile(x, 100. * q)


#print d.shape

t0 = 2456836.2964455
per =  0.925545613
t = np.linspace(t0, t0 + 1.0*per, 1000)
rp = 0.1146
stellar_grid = spiderman.stellar_grid.gen_grid(1.1e-6, 1.7e-6, stellar_model = "blackbody")
dilution = 0.12777
Ts = 6110.

phs = (t - t0)/per
n = 100
peak_phases = np.zeros(n)

#sph = [2.32996923e+03, 1.61955101e+02, 0.00000000e+00, 5.90872619e+02]
#sph = [2232.382577585771, 16.68482360422172, 1000.0, 580.3525685355843]
sph = [2114., 16.3, 531., 640.]
lc = spiderman_spherical.lc(t, rp, Ts, 1.1e-6, 1.7e-6, sph, dilution = dilution, eclipse = False, stellar_grid = stellar_grid)

ind = lc == lc.max()
print "peak phase = ", phs[ind]

plt.axvline(0.5)
plt.plot(phs, lc)
plt.show()


#for i in range(len(t)): print phs[i], (lc[i] - 1.)*1.e3


