import numpy as np

#made the file teffs.txt based on output from phased_espec.py
teffs = np.genfromtxt("teffs.txt")

print "\\begin{deluxetable}{LLL}"
print "\\tablecolumns{3}"
print "\\tablewidth{0pt}:"
print "\\tablecaption{Phase-resolved Effective Temperature \label{table:teffs}}"
print "\\tablehead{"
print "\colhead{Orbital Phase} & \colhead{$T_\mathrm{eff}$} & \colhead{$\chi^2_\\nu$}}"
print "\startdata"

T = teffs[:,0]
err = teffs[:,1]
chi2 = teffs[:,2]
phase = teffs[:,3]
phase1 = np.array([0.06, 0.15, 0.25, 0.35, 0.44, 0.56, 0.65, 0.75, 0.85])
phase2 = np.array([0.15, 0.25, 0.35, 0.44, 0.56, 0.65, 0.75, 0.85, 0.94])

for i in range(len(teffs)): print phase1[i], "-", phase2[i],  "&", int(T[i]), "\\pm", int(err[i]), "&", chi2[i],  "\\\\"

print "\enddata"
print "\\vspace{-0.8cm}"
print "\\tablecomments{comments}"
print "\end{deluxetable}"

