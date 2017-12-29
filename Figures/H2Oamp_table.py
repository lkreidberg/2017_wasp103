import numpy as np

h2o_amp = np.genfromtxt("h2o_amp_1.15_1.3_1.35_1.55.txt")
amp = h2o_amp[:,0]
amp_e = h2o_amp[:,1]

print "\\begin{deluxetable}{llllll}"
print "\\tablecolumns{6}"
print "\\tablewidth{0pt}:"
print "\\tablecaption{Water feature amplitude \label{table:h2o_amp}}"
print "\\tablehead{"
print "\colhead{\,} & \colhead{Object} & \colhead{$T_\mathrm{eff}$} & \colhead{logg} & \colhead{Wide Amp} & \colhead{Narrow Amp}}"
print "\startdata"

names = ["W103b night", "W103b quadrature", "W103b dayside", "1320+0409", "0428-2253", "0003-2822", "1RXSJ160929.1", "J1610-1913B", "TWA22A"]
objects = ["Hot Jupiter", "\,", "\,", "Brown Dwarf", "\,", "\,", "Imaged Companion", "\,", "\,"]
Teff = [1922, 2408, 2982, 1875, 2429, 2889, 1800, 2400, 3000]
Teff_e = [0, 0, 0, 70, 80, 80, 0, 0, 0] 
logg = [3.20, 3.20, 3.20, 5.19, 5.22, 5.18, 0, 0, 0]
logg_e = [0.04, 0.04, 0.04, 0.16, 0.09, 0.04, 0, 0, 0]

#category (e.g. "BDs"), object, H2O amp1, H2O amp2, teff, logg, rotation period?

n = 9           #nine objects

for i in range(n): print objects[i], "&", names[i], "&", "$", int(Teff[i]), "\\pm", int(Teff_e[i]), "$", "&", "$", logg[i], "\\pm", logg_e[i], "$", "&", "$", amp[i], "\\pm", amp_e[i], "$", "&", "$", amp[i], "\\pm", amp_e[i], "$",  "\\\\"

print "\enddata"
print "\\vspace{-0.8cm}"
print "\\tablecomments{comments}"
print "\end{deluxetable}"

