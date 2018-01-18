import numpy as np

wfc3 = np.genfromtxt("teffs.txt")
spt = np.genfromtxt("teffs_ch2.txt")


for i in range(len(wfc3)): print (wfc3[i,0] - spt[i,0])/np.sqrt(wfc3[i,1]**2 + spt[i,1]**2)
