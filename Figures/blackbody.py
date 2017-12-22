import numpy as np
import matplotlib.pyplot as plt
from astropy import modeling

w = np.linspace(1, 2, 10)
flux = modeling.blackbody.blackbody_lambda(w*1.e4, 6100.)       #wavelength in angstroms
# return units: math:`erg \\; cm^{-2} s^{-1} \\mathring{A}^{-1} sr^{-1}

flux = flux*4.*np.pi*1.e4*w        #multiplies by 4 pi steradians * 10000. angstroms/micron * micron

scale = 1.863e17  #(470 pc/jupiter radius)^2
plt.plot(w, flux/scale*1.e-3)

plt.show()
