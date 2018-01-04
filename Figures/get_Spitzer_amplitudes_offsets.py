import numpy as np
import pickle
import matplotlib.pyplot as plt

#return 1.+p.amp1[v_num]*np.cos(2.*np.pi*(t-p.theta1[v_num])/p.per[v_num])  //HST

#IND(y,i) = c1a*cos(2*pi*(IND(x,i)-c1o)/p) + c                              //Spitzer

#p = 0.925545613
p = 9.2554e-01
depth, depth_err =  4.85018804e-03, 2.16170975e-04
c1a, c1a_err =  2.00086275e-03, 9.82717991e-05
c1o, c1o_err =  5.32560819e-01, 6.58478628e-03

toffset = 2457162.

offset = (((c1o + toffset)/p)%1)*180./np.pi
phase_offset = ((c1o+toffset)/p)%1

#print "CH2 offset in degrees", offset, "+/-", c1o_err*180./np.pi
#print "Ch2 offset in phase", phase_offset
#print "is there truncation error in the period? this is a problem for toffset >> p"

pic = pickle.load(open("spitzer_ch2_sine_bestfit.pic", "rb"))

plt.subplot(211)
plt.plot(pic[0], pic[3])
x = np.linspace(0, p, 100)
#plt.plot(x/p, 1. - c1a*np.cos(2*np.pi*(x-phase_offset*p)/p)) 

#guess is more accurate!!
guess2 = 0.017
y =  1. - c1a*np.cos(2*np.pi*(x+guess2)/p)
plt.plot(x/p, 1. - c1a*np.cos(2*np.pi*(x+guess2)/p)) 
print "ch2 offset by eye = ", guess2/p*180./np.pi, "+/-", c1o_err*180./np.pi
print "ch2 amplitude = ", c1a*2./depth, "+/-", c1a*2./depth*np.sqrt((c1a_err/c1a)**2 + (depth_err/depth)**2) 
#print "ch2 guess amplitude = ", (y.max() - y.min())/depth  #it's exactly the same

plt.ylim(0.997, 1.003)


#Ch 1
depth, depth_err =  3.68401784e-03, 2.99006690e-04
c1a, c1a_err = 1.84589958e-03, 2.37261751e-04
c1o, c1o_err = 1.46311670e-01, 1.17483920e-02

toffset = 2457170.7
offset = (((c1o + toffset)/p)%1)*180./np.pi
phase_offset = ((c1o+toffset)/p)%1
#print "CH1 offset in degrees", offset 
#print "CH1 offset in phase", phase_offset

plt.subplot(212)
pic = pickle.load(open("Ch1_best_fits/2017-09-15_11:33-sincos2/bestfit.pic", "rb"))
plt.plot(pic[0], pic[3])
x = np.linspace(0, p, 100)
guess1 = 0.032
plt.plot(x/p, 1. - c1a*np.cos(2*np.pi*(x+guess1)/p)/1.17)       #divide by 1.17 for dilution correction (this was run with dilution included??)
plt.ylim(0.997, 1.003)

print "ch1 offset by eye = ", guess1/p*180./np.pi, "+/-", c1o_err*180./np.pi
print "ch1 amplitude = ", c1a*2./depth/1.17, "+/-", c1a*2./1.17/depth*np.sqrt((c1a_err/c1a)**2 + (depth_err/depth)**2) 
 
plt.show()
