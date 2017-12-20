import numpy as np
import matplotlib.pyplot as plt

def blackbody(l,T): 
	h = 6.62607e-34   #J/s
        c = 2.997925e8       #m/s
        k = 1.38065e-23     #J/K
        return 2*h*c**2/(l**5*(np.exp(h*c/(l*k*T))-1))                 #[W/sr/m^3]

def best_fit_bb(w, y, e, rprs):
	Ts = np.linspace(300, 3600, 300)
	chibest = 10000.
	Tbest = 0.	
	Tstar = 6110.
	w = np.array(w)
	waves_hires = np.linspace(1.0, 2.0, 100)


	#get stellar spectrum
	star = np.genfromtxt("wasp103_sed_fluxes.out")
	star_bb = np.interp(w, star[:,0], star[:,1])*1.e24/(w*np.pi*4.)
        #Tstar = 6110.
        #star_bb = blackbody(w*1.0e-6, Tstar)
	outliers = 0.

	for T in Ts:
		model = blackbody(w*1.0e-6, T)/star_bb*rprs**2
		chi2 = np.sum((y - model)**2/e**2)
		if chi2 < chibest: 
			chibest, Tbest, outliers = chi2, T, (y-model)/e
	star_bb_hires = np.interp(waves_hires, star[:,0], star[:,1])*1.e24/(waves_hires*np.pi*4.)
        #star_bb_hires = blackbody(waves_hires*1.0e-6, Tstar)
	#print "Best fit blackbody temp, chisq, and outliers: ", Tbest, chibest/(len(e)-1), outliers
	print "Best fit blackbody temp, chi2: ", Tbest, chibest/(len(y) - 1.)
	return waves_hires, blackbody(waves_hires*1.0e-6, Tbest)/star_bb_hires*rprs**2

d = np.genfromtxt("temp_espec.txt")
waves = d[:,0]
fpfs, fp_err = d[:,1], d[:,2]

rprs = 0.1127
wave_hires, model_hires = best_fit_bb(waves, fpfs, fp_err, rprs)

plt.plot(wave_hires, model_hires, color='0.5', label='blackbody', linestyle = 'dashed', zorder = -20)
plt.errorbar(waves, fpfs, fp_err, fmt = '.k')
plt.xlim(1,1.8)
plt.ylim(1e-3, 1.8e-3)
plt.xlabel("wavelength (microns")
plt.ylabel("fp/fs")
plt.show()
