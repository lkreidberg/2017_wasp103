import pickle
import numpy as np

z = pickle.load(open("zhang.p", "rb"))
s = pickle.load(open("spherical.p", "rb"))
h = pickle.load(open("hotspot_t.p", "rb"))
sc = pickle.load(open("sincos.p", "rb"))

s, se  = s[1][1::], s[2][1::]
z, ze  = z[1][1::], z[2][1::]
h, he  = h[1][1::], h[2][1::]
sc, sce  = sc[1][1::], sc[2][1::]

zdiff = np.abs((s - z)/np.sqrt(se**2 + ze**2))
hdiff = np.abs((s - h)/np.sqrt(he**2 + se**2))
sdiff = np.abs((z - h)/np.sqrt(he**2 + ze**2))

print "median deviation", "fraction of bins within one sigma"
print "zhang - spherical", np.median(zdiff), 1.*sum(zdiff < 1.)/len(zdiff)
print "spherical - hotspot", np.median(hdiff), 1.*sum(hdiff < 1.)/len(hdiff)
print "zhang - hotspot", np.median(sdiff), 1.*sum(sdiff < 1.)/len(sdiff)


diff = np.abs((z - sc)/np.sqrt(ze**2 + sce**2))
print "zhang - sincos", np.median(diff), 1.*sum(diff < 1.)/len(diff)

diff = np.abs((s - sc)/np.sqrt(se**2 + sce**2))
print "spherical - sincos", np.median(diff), 1.*sum(diff < 1.)/len(diff)

diff = np.abs((h - sc)/np.sqrt(he**2 + sce**2))
print "hotspot - sincos", np.median(diff), 1.*sum(diff < 1.)/len(diff)



print "hostpot  - sincos all", diff, np.mean(diff[24:48]), len(diff)
