import numpy as np

h2o_amp1 = np.genfromtxt("h2o_amp_1.15_1.3_1.35_1.5.txt")
h2o_amp2 = np.genfromtxt('h2o_amp_1.2_1.35_1.5_1.65.txt')
amp1 = h2o_amp1[:,0]
amp1_e = h2o_amp1[:,1]
amp2 = h2o_amp2[:,0]
amp2_e = h2o_amp2[:,1]


names = ["W103b night", "W103b quadrature", "W103b dayside", "2MASS J1320+0409", "2MASS J0428-2253", "2MASS J0003-2822", "CD-35 2722\\tablenotemark{3}", "USco 1610-1913B\\tablenotemark{2}", "TWA 22A\\tablenotemark{3}"]


objects = ["Hot Jupiter", "\,", "\,", "Brown Dwarf", "\,", "\,", "Imaged Companion", "\,", "\,"]
Teff = [1922, 2408, 2982, 1875, 2429, 2889, 1800, 2400, 3000]
Teff_e = [40, 20, 10, 70, 80, 80, 100, 150, 100.] 
logg = [3.20, 3.20, 3.20, 5.19, 5.22, 5.18, 4.5, 0, 4.5]
logg_e = [0.04, 0.04, 0.04, 0.16, 0.09, 0.04, 0.5, 0, 0.5]

##1RXS: Lafreniere 2010
#CD-35: Wahhaj 2011
#USco 1610: Aller 2013
#TWA 22A: Bonnefoy et al. 2014 Table 6

#category (e.g. "BDs"), object, H2O amp1, H2O amp2, teff, logg, rotation period?

n = 9           #nine objects

print "\\begin{deluxetable*}{llCCCC}"
print "\\tablecolumns{6}"
print "\\tablewidth{0pt}:"
print "\\tablecaption{Source Properties \label{table:sources}}"
print "\\tablehead{"
print "\colhead{\,} & \colhead{Object} & \colhead{$T_\mathrm{eff}$} & \colhead{logg} & \colhead{H$_2$O A$_1$} & \colhead{H$_2$O A$_2$}}"
print "\startdata"
#for i in range(n): print objects[i], "&", names[i], "&", "$", int(Teff[i]), "\\pm", int(Teff_e[i]), "$", "&", "$", logg[i], "\\pm", logg_e[i], "$", "&", "$", amp1[i], "\\pm", amp1_e[i], "$", "&", "$", amp2[i], "\\pm", amp2_e[i], "$",  "\\\\"
for i in range(n): print objects[i], "&", names[i], "&",  int(Teff[i]), "\\pm", int(Teff_e[i]), "&",  logg[i], "\\pm", logg_e[i],  "&",  '{0:0.2f}'.format(amp1[i]), "\\pm", '{0:0.1e}'.format(amp1_e[i]),  "&", '{0:0.2f}'.format(amp2[i]), "\\pm", '{0:0.1e}'.format(amp2_e[i]),  "\\\\"

print "\enddata"
print "\end{deluxetable*}"

