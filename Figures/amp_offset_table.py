import numpy as np


print "\\begin{deluxetable}{LLL}"
print "\\tablecolumns{3}"
print "\\tablewidth{0pt}:"
print "\\tablecaption{Phase Curve Properties \label{table:amps_offsets}}"
print "\\tablehead{"
print "\colhead{Wavelength} & \colhead{Amplitude} & \colhead{Offset}}"
print "\startdata"


d = np.genfromtxt("WFC3_amplitudes_offsets.txt")
#columns: 1) phase shift, 2-3)phase shift err, 4) amp, 5-6) amp err, 7) fpfs, 8-9) fpfs_err

off, off_err = d[:,1]*180./np.pi, (d[:,2]+ d[:,3])/2.*180./np.pi
amp, amp_err = d[:,4], (d[:,5]+d[:,6])/2.
fpfs, fpfs_err = d[:,7], (d[:,8]+d[:,9])/2.
wave = np.linspace(1.175, 1.625, 10)

for i in range(len(d)): print wave[i], "&",  '{0:0.2f}'.format(amp[i]), "\\pm", '{0:0.2f}'.format(amp_err[i]), "&", '{0:0.1f}'.format(off[i]), "\\pm", '{0:0.1f}'.format(off_err[i]),  "\\\\"

print "\enddata"
print "\end{deluxetable}"

