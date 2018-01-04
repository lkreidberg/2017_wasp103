import numpy as np

def weighted_mean(data, err):                #calculates the weighted mean for data points data with variances err
    weights = 1.0/err
    mu = np.sum(data*weights)/np.sum(weights)
    var = 1.0/np.sum(weights)
    return [mu, var]                #returns weighted mean and variance

d = np.genfromtxt("WFC3_amplitudes_offsets.txt")
#columns: 1) phase shift, 2-3)phase shift err, 4) amp, 5-6) amp err, 7) fpfs, 8-9) fpfs_err

off, off_err = -1.*d[:,1]*180./np.pi, (d[:,2]+ d[:,3])/2.*180./np.pi
amp, amp_err = d[:,4], (d[:,5]+d[:,6])/2.

fpfs, fpfs_err = d[:,7], (d[:,8]+d[:,9])/2.
wave = np.linspace(1.175, 1.625, 10)

amp, amp_err = weighted_mean(amp, amp_err) 
off, off_err = weighted_mean(off, off_err)

gcm_wfc3 = np.genfromtxt("gcm_amp_off_wfc3.txt")
models = ["nominal GCM", "[Fe/H] = 0.5 GCM", "$\\tau_\mathrm{drag4}$ GCM", "$\\tau_\mathrm{drag3}$ GCM"]

print "\\begin{deluxetable}{llLL}"
print "\\tablecolumns{4}"
print "\\tablewidth{0pt}:"
print "\\tablecaption{Phase Curve Properties \label{table:amps_offsets}}"
print "\\tablehead{"
print "\colhead{Bandpass} & \colhead{Source} & \colhead{Amplitude} & \colhead{Offset)}\\\\"
print "\colhead{\,} & \colhead{\,} & \colhead{\,} & \colhead{(Degrees)}}"
print "\startdata"

#WFC3
print "WFC3", "&", "data", "&",  '{0:0.2f}'.format(amp), "\\pm", '{0:0.2f}'.format(amp_err), "&", '{0:0.1f}'.format(off), "\\pm", '{0:0.1f}'.format(off_err),  "\\\\"

for i, model in enumerate(models): print "\,", "&", model, "&",  '{0:0.2f}'.format(gcm_wfc3[i, 2]), "&", '{0:0.2f}'.format(gcm_wfc3[i,1]),   "\\\\"


#Spitzer data
wave = np.array([3.6, 4.5])
off = np.array([2.0, 1.0])
off_err = np.array([0.7, 0.4])
amp = np.array([0.86, 0.83])
amp_err = np.array([0.13, 0.05])

gcm_ch1 = np.genfromtxt("gcm_amp_off_ch1.txt")
gcm_ch2 = np.genfromtxt("gcm_amp_off_ch2.txt")

#Spitzer Ch 1
i = 0
print "Spitzer 3.6 $\mu$m", "&", "data", "&",  '{0:0.2f}'.format(amp[i]), "\\pm", '{0:0.2f}'.format(amp_err[i]), "&", '{0:0.1f}'.format(off[i]), "\\pm", '{0:0.1f}'.format(off_err[i]),  "\\\\"

for i, model in enumerate(models): print "\,", "&", model, "&",  '{0:0.2f}'.format(gcm_ch1[i, 2]), "&", '{0:0.2f}'.format(gcm_ch1[i,1]),   "\\\\"


#Spitzer Ch 2
i = 1
print "Spitzer 4.5 $\mu$m", "&", "data", "&",  '{0:0.2f}'.format(amp[i]), "\\pm", '{0:0.2f}'.format(amp_err[i]), "&", '{0:0.1f}'.format(off[i]), "\\pm", '{0:0.1f}'.format(off_err[i]),  "\\\\"

for i, model in enumerate(models): print "\,", "&", model, "&",  '{0:0.2f}'.format(gcm_ch2[i, 2]), "&", '{0:0.2f}'.format(gcm_ch2[i,1]),   "\\\\"

print "\enddata"
print "\end{deluxetable}"

