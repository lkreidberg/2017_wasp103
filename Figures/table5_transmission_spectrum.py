import numpy as np

#made the file teffs.txt based on output from phased_espec.py
t0 = np.genfromtxt("tspec_tnight0.txt")
t1700 = np.genfromtxt("tspec_tnight1700.txt")


print "\\begin{deluxetable}{LLLL}"
print "\\tablecolumns{4}"
print "\\tablewidth{0pt}"
print "\\tablecaption{WASP-103b Transmission Spectrum\label{table:tspec}}"
print "\\tablehead{"
print "\colhead{Wavelength} & \colhead{$(R_p/R_s)^2$ (\%)} & \colhead{$(R_p/R_s)^2$ (\%)} & \colhead{Error}\\\\"
print "\colhead{(micron)} & \colhead{($T_\mathrm{night} = 0$ K)} & \colhead{($T_\mathrm{night} = 1700$ K)} & \colhead{(\%)}}"
print "\startdata"

#for i in range(len(t0)): print t0[i,0],  "&", '{0:0.4f}'.format(100.*t0[i,1]), "\\pm", '{0:0.4f}'.format(100.*t0[i,2]), "&", '{0:0.4f}'.format(t1700[i,1]*100.), "\\pm", '{0:0.4f}'.format(t1700[i,2]*100.),  "\\\\"
for i in range(len(t0)): print t0[i,0],  "&", '{0:0.4f}'.format(100.*t0[i,1]), "&", '{0:0.4f}'.format(t1700[i,1]*100.), "&", '{0:0.4f}'.format(t1700[i,2]*100.),  "\\\\"

print "\enddata"
print "\\vspace{-0.8cm}"
print "\\tablecomments{comments}"
print "\end{deluxetable}"

