import numpy as np

def weighted_mean(data, err):                #calculates the weighted mean for data points data with variances err
    weights = 1.0/err
    mu = np.sum(data*weights)/np.sum(weights)
    var = 1.0/np.sum(weights)
    return [mu, var]                #returns weighted mean and variance


"\\begin{deluxetable}{LLL}"
print "\\tablecolumns{3}"
print "\\tablewidth{0pt}:"
print "\\tablecaption{Phase Curve Properties \label{table:amps_offsets}}"
print "\\tablehead{"
print "\colhead{Wavelength} & \colhead{Amplitude} & \colhead{Offset (Degrees)}}"
print "\startdata"


d = np.genfromtxt("WFC3_amplitudes_offsets.txt")
#columns: 1) phase shift, 2-3)phase shift err, 4) amp, 5-6) amp err, 7) fpfs, 8-9) fpfs_err

off, off_err = -1.*d[:,1]*180./np.pi, (d[:,2]+ d[:,3])/2.*180./np.pi
amp, amp_err = d[:,4], (d[:,5]+d[:,6])/2.

fpfs, fpfs_err = d[:,7], (d[:,8]+d[:,9])/2.
wave = np.linspace(1.175, 1.625, 10)

print "WFC3 average amplitude = ", weighted_mean(amp, amp_err) 
print "WFC3 average offset = ", weighted_mean(off, off_err)

for i in range(len(d)): print wave[i], "&",  '{0:0.2f}'.format(amp[i]), "\\pm", '{0:0.2f}'.format(amp_err[i]), "&", '{0:0.1f}'.format(off[i]), "\\pm", '{0:0.1f}'.format(off_err[i]),  "\\\\"

#Spitzer data
wave = np.array([3.6, 4.5])
off = np.array([2.0, 1.0])
off_err = np.array([0.7, 0.4])
amp = np.array([0.86, 0.83])
amp_err = np.array([0.13, 0.05])

for i in range(2): print wave[i], "&",  '{0:0.2f}'.format(amp[i]), "\\pm", '{0:0.2f}'.format(amp_err[i]), "&", '{0:0.1f}'.format(off[i]), "\\pm", '{0:0.1f}'.format(off_err[i]),  "\\\\"

print "\enddata"
print "\end{deluxetable}"

