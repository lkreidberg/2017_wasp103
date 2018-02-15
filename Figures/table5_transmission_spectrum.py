import numpy as np

#made the file teffs.txt based on output from phased_espec.py
t0 = np.genfromtxt("tspec_tnight0.txt")
t1700 = np.genfromtxt("tspec_tnight1700.txt")


print "\\begin{deluxetable}{LLL}"
print "\\tablecolumns{3}"
print "\\tablewidth{0pt}"
print "\\tablecaption{WASP-103b Transmission Spectrum\label{table:tspec}}"
print "\\tablehead{"
print "\colhead{Wavelength} & \colhead{Transit Depth (\%)} & \colhead{Transit Depth (\%)} \\\\"
print "\colhead{(micron)} & \colhead{($T_\mathrm{night} = 0$ K)} & \colhead{($T_\mathrm{night} = 1700$ K)}}"
print "\startdata"

for i in range(len(t0)): print t0[i,0],  "&", '{0:0.6f}'.format(100.*t0[i,1]), "\\pm", '{0:0.6f}'.format(t0[i,2]), "&", '{0:0.6f}'.format(t1700[i,1]), "\\pm", '{0:0.6f}'.format(t1700[i,2]),  "\\\\"

print "\enddata"
print "\\vspace{-0.8cm}"
print "\\tablecomments{comments}"
print "\end{deluxetable}"

