import matplotlib.pyplot as plt
import numpy as np
import os, glob

path1 = "phys_spectra"
path2 = "sine_spectra"
path3 = "phys_spectra_nov2"

files1 = glob.glob(os.path.join(path1, "*.txt"))		
files2 = glob.glob(os.path.join(path2, "*.txt"))		
files3 = glob.glob(os.path.join(path3, "*.txt"))		

deviate = []

max = 6
for i, f in enumerate(files1):
	plt.clf()
	p  = np.genfromtxt(f)
	s = np.genfromtxt(files2[i])
	x = np.genfromtxt(files3[i])
	plt.errorbar(p[0:max,0], p[0:max,1], p[0:max,2], fmt = '.k')
	plt.errorbar(s[0:max,0], s[0:max,1], s[0:max,2], fmt = '.r')
	plt.errorbar(x[0:max,0], x[0:max,1], x[0:max,2], fmt = '.g')
	for j in range(len(p[0:max,1])-2): deviate.append((p[j,1] - s[j,1])/s[j,2])
	plt.xlim(1.1, 1.7)
	figname = "compare_" + "{0:0.1f}".format(p[0,4] - 0.05) + ".png"	
	plt.savefig(figname)


deviate = np.array(deviate)
print "mean and median of deviates", np.mean(deviate), np.median(deviate)
