import numpy as np
import pickle
import matplotlib.pyplot as plt
from lc_fit import Model, LightCurveData, PrintParams
import spiderman

p = pickle.load(open("WFC3_1.18_bestfit.pic", "rb"))
per = 0.925545613

d, m, par = p[0], p[1], p[2]		
ind = d.err < 9.0e7			#indices of outliers
    
err = d.err[ind]/d.flux[ind]		#normalized per point uncertainty 
phase = m.phase[ind]
data_corr = m.data_corr[ind]
bestfit = m.bestfit[ind]

ind = np.argsort(phase)
phase, bestfit = phase[ind], bestfit[ind]
plt.plot(phase, bestfit)

amp1 =  par[d.par_order['amp1']*d.nvisit:(1 + d.par_order['amp1'])*d.nvisit][0]
theta1 = par[d.par_order['theta1']*d.nvisit:(1 + d.par_order['theta1'])*d.nvisit][0]
depth = par[d.par_order['fp']*d.nvisit:(1 + d.par_order['fp'])*d.nvisit][0]
print theta1/per*180./np.pi

theta1, amp1, depth  = -0.0313997651933, 1.1066970271, 0.000519807824444
print theta1/per*180./np.pi
#print 'offset', theta1/per*180./np.pi
#print 'offset2', theta1*180./np.pi

t = np.linspace(0, per, 100)
plt.plot(t/per, 1.+ depth*(1. - amp1*np.cos(2.*np.pi*(t-theta1)/per)))

plt.ylim(0.999, 1.0012)
plt.show()
