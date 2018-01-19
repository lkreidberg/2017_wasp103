import pickle
import numpy as np

z = pickle.load(open("zhang.p", "rb"))
s = pickle.load(open("spherical.p", "rb"))
h = pickle.load(open("hotspot_t.p", "rb"))

s, se  = s[1][1::], s[2][1::]
z, ze  = z[1][1::], z[2][1::]
h, he  = h[1][1::], h[2][1::]

zdiff = np.abs((s - z)/np.sqrt(s**2 + z**2))
hdiff = np.abs((s - h)/np.sqrt(h**2 + s**2))
sdiff = np.abs((z - h)/np.sqrt(h**2 + z**2))

print np.median(zdiff), 1.*sum(zdiff < 1.)/len(zdiff)
print np.median(hdiff), 1.*sum(hdiff < 1.)/len(hdiff)
print np.median(sdiff), 1.*sum(sdiff < 1.)/len(sdiff)
