import corner
import numpy as np
import matplotlib.pyplot as plt
import os, glob

path = "../ECLIPSE_MCMC/"
files = glob.glob(os.path.join(path, "mcmc_out*npy"))		

for f in files:
	figname = "pairs"+f[25:29]+".png"
	print figname
	
	d = np.load(f)
	labels = ["fp", "c", "c", "c", "c", "v", "v", "v", "v", "r1", "r2", "scale"]
	fig = corner.corner(d, labels = labels, show_titles=True)
	plt.savefig(figname)	
