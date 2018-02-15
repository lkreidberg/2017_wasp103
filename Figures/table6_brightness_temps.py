import numpy as np

#made the file teffs.txt based on output from phased_espec.py
teffs = np.genfromtxt("teffs.txt")
teffs_s1 = np.genfromtxt("teffs_ch1.txt")
teffs_s2 = np.genfromtxt("teffs_ch2.txt")

print "\\begin{deluxetable}{LLLLL}"
print "\\tablecolumns{5}"
print "\\tablewidth{0pt}"
print "\\tablecaption{Phase-resolved Brightness Temperatures \label{table:teffs}}"
print "\\tablehead{"
print "\colhead{Orbital Phase} & \colhead{$T_\mathrm{b}$} & \colhead{$T_\mathrm{b}$} & \colhead{$T_\mathrm{b}$} & \colhead{$\chi^2_\\nu$} \\\\"
print "\colhead{\,} & \colhead{WFC3} & \colhead{Ch. 1} & \colhead{Ch. 2} & \colhead{(9 dof)}}"
print "\startdata"

T = teffs[:,0]
err = teffs[:,1]
chi2 = teffs[:,2]
phase = teffs[:,3]
phase1 = np.array([0.06, 0.15, 0.25, 0.35, 0.44, 0.56, 0.65, 0.75, 0.85])
phase2 = np.array([0.15, 0.25, 0.35, 0.44, 0.56, 0.65, 0.75, 0.85, 0.94])

T1 = teffs_s1[:,0]
err1 = teffs_s1[:,1]
T2 = teffs_s2[:,0]
err2 = teffs_s2[:,1]

for i in range(len(teffs)): print phase1[i], "-", phase2[i],  "&", int(T[i]), "\\pm", int(err[i]), "&", int(T1[i]), "\\pm", int(err1[i]), "&", int(T2[i]), "\\pm", int(err2[i]), "&", chi2[i],  "\\\\"

print "\enddata"
print "\\vspace{-0.8cm}"
print "\\tablecomments{comments}"
print "\end{deluxetable}"

