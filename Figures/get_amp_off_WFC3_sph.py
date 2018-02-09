import numpy as np
import spiderman_spherical
import spiderman

d = np.load("/Users/lkreidberg/Desktop/Projects/Observations/HST/WASP103_HST_all/SPHERICAL_MCMC_020818/flatchain2_2018_02_08_16:36:5766532957.2389.npy")

#print d.shape

t = np.linspace(0, 0.92554, 100)
rp = 0.1146
stellar_grid = spiderman.stellar_grid.gen_grid(1.1e-6, 1.7e-6, stellar_model = "/Users/lkreidberg/Desktop/Documents/Papers/2017_wasp103/Figures/w103_spectrum.txt")
dilution = 0.12777
Ts = 6110.

for i in range(1):
	sph = [d[i, 16], d[i,17], 0., d[i,18]]
	lc = spiderman_spherical.lc(t, rp, Ts, 1.1e-6, 1.7e-6, sph, dilution = dilution, eclipse = False, stellar_grid = stellar_grid)
	ind = lc == lc.max()
	#print t[ind]

for i in range(len(t)): print t[i], lc[i]
