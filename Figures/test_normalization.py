import numpy as np
import matplotlib.pyplot as plt


path = "BD_spectra/"
files = ["0024-0158 (M9.5) SED.txt", "0428-2253 (L0.5) SED.txt"]

distance = np.array([11.55, 25.99])
label = ["0024-0158", "0428-2253"]

for i, f in enumerate(files):

	d = np.genfromtxt(path + f)
	path = "BD_spectra/"

	plt.plot(d[:,0], d[:,1]*distance[i]**2, label = label[i])
	
plt.legend()
plt.ylabel("flux density $\\times$ distance$^2$")
plt.xlabel("wavelength (micron)")
plt.xlim(0, 3)

plt.show()


