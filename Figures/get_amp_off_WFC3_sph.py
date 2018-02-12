import numpy as np
import spiderman_spherical
import spiderman

def quantile(x, q): return np.percentile(x, 100. * q)

amp = False
offset = True

d = np.load("/Users/lkreidberg/Desktop/Projects/Observations/HST/WASP103_HST_all/SPHERICAL_MCMC_020818/flatchain2_2018_02_08_16:36:5766532957.2389.npy")

#print d.shape

t0 = 2456836.2964455
per =  0.925545613
t = np.linspace(t0, t0 + 1.0*per, 1000)
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
phs = (t - t0)/per
n = 100
peak_phases, amplitudes = np.zeros(n), np.zeros(n)

print "min/max sph0", np.min(d[:,16]), np.max(d[:,16])

for i in range(n):
	sph = [d[i, 16], d[i,17], 0., d[i,18]]
	print sph
	if offset:
		lc = spiderman_spherical.lc(t, rp, Ts, 1.1e-6, 1.7e-6, sph, dilution = dilution, eclipse = False, stellar_grid = stellar_grid)
		ind = lc == lc.max()
		peak_phases[i] = phs[ind]
	if amp:
		lc = spiderman_spherical.lc(t, rp, Ts, 1.1e-6, 1.7e-6, sph, dilution = dilution, eclipse = True, stellar_grid = stellar_grid)
		ind = lc == lc.max()
		peak_phases[i] = phs[ind]

print quantile(peak_phases, np.array([0.16, 0.5, 0.84]))

#for i in range(len(t)): print phs[i], (lc[i] - 1.)*1.e3

