import numpy as np

teffs = np.genfromtxt("teffs.txt")

print "\\begin{deluxetable}{lll}"
print "\\tablecolumns{3}"
print "\\tablewidth{0pt}:"
print "\\tablecaption{Phase-resolved Effective Temperature \label{table:teffs}}"
print "\\tablehead{"
print "\colhead{Phase} & \colhead{Effective Temperature} & \colhead{$\chi^2_\\nu$}}"
print "\startdata"

T = teffs[:,0]
err = teffs[:,1]
chi2 = teffs[:,2]
phase = teffs[:,3]

for i in range(len(teffs)): print phase[i], "&", "$", int(T[i]), "\\pm", int(err[i]), "$", "&", chi2[i],  "\\\\"

print "\enddata"
print "\\vspace{-0.8cm}"
print "\\tablecomments{comments}"
print "\end{deluxetable}"

