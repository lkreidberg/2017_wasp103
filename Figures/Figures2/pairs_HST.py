import corner
import numpy as np
import matplotlib.pyplot as plt
import os, glob

path = "../PHYSICAL_MODEL_MCMC/"
files = glob.glob(os.path.join(path, "mcmc_out*npy"))		

for f in files:
	figname = "pairs"+f[32:36]+".png"
	print figname
	
	d = np.load(f)
	labels = ["rp", "u1", "c", "c", "c", "c", "v", "v", "v", "v", "v2", "v2", "r1", "r2", "scale", "xi", "delta_T", "T_n", "T_s", "x", "x", "x", "x"]
	fig = corner.corner(d, labels = labels, show_titles=True)
	plt.savefig(figname)	
