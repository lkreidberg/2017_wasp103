import numpy as np
import spiderman_spherical
import spiderman

def quantile(x, q): return np.percentile(x, 100. * q)

d = np.load("/Users/lkreidberg/Desktop/Projects/Observations/Spitzer/WASP103_11099/Ch2/analysis/run/fgc/ap2750715/fit_phase_curve/2018-02-08_16:51-spherical_longmcmc/d-WA103bo21-allparams-trq1lnspiderman_sphericalbli.npy") 

d = d.T
print d[:,0].shape

t0 = 2456836.2964455
per =  0.925545613
t = np.linspace(t0, t0 + 1.0*per, 1000)
rp = 0.1146
stellar_grid = spiderman.stellar_grid.gen_grid(1.1e-6, 1.7e-6, stellar_model = "/Users/lkreidberg/Desktop/Documents/Papers/2017_wasp103/Figures/w103_spectrum.txt")
dilution = 0.12777
Ts = 6110.

phs = (t - t0)/per
n = 1000
peak_phases = np.zeros(n)

print "min/max sph0", np.min(d[:,27]), np.max(d[:,27])

for i in range(n):
	sph = [d[i, 27], d[i,28], 0., d[i,30]]
	lc = spiderman_spherical.lc(t, rp, Ts, 1.1e-6, 1.7e-6, sph, dilution = dilution, eclipse = False, stellar_grid = stellar_grid)
	ind = lc == lc.max()
	peak_phases[i] = phs[ind]

qs = quantile(peak_phases, np.array([0.16, 0.5, 0.84]))
print "max phase", np.max(peak_phases)
print (qs[1] - 0.5)*180/np.pi, (qs[1] - qs[0])*180./np.pi, (qs[2] - qs[1])*180./np.pi

#for i in range(len(t)): print phs[i], (lc[i] - 1.)*1.e3

