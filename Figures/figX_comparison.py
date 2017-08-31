import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astropy.modeling.blackbody import blackbody_lambda
import pyfits
from astropy import units as u

plt.figure()
gs = gridspec.GridSpec(4, 3, hspace=0.05, wspace=0.25)

#Brown dwarf spectra
"""temp object
1918 +- 72  2MASSJ13204427p0409045_spec_app.txt
1988 +- 75  2MASSIJ2104149m103736_spec_app.txt
1995 +- 77  2MASSWJ1506544p132106_spec_app.txt
2212 +- 55  2MASPJ0345432p254023_spec_app.txt
2275 +- 56  LHS2924_spec_app.txt
2585 +- 46  2MASSJ03205965p1854233_spec_app.txt
2658 +- 121 LHS2021_spec_app.txt
2817 +- 59  2MASSJ03510004m0052452_spec_app.txt
2873 +- 75  2MASSJ00034227m2822410_spec_app.txt"""

path = "BD_spectra/"
files = ["2MASSJ13204427p0409045_spec_app.txt", "2MASPJ0345432p254023_spec_app.txt", "2MASSJ03205965p1854233_spec_app.txt", "2MASSJ00034227m2822410_spec_app.txt"]

for i, f in enumerate(files):

    ax = plt.subplot(gs[i,0])

    wavelengths = [1, 2] * u.micron
    temperature = 2000 * u.K
    flux = blackbody_lambda(wavelengths, temperature)
    plt.plot(wavelengths, flux)

    d = np.genfromtxt(path + f)
    plt.plot(d[:,0], d[:,1])
    plt.xlim(1,2)


"""#Directly imaged spectra
2M1207 A     2500  100 3.5  0.5 2640
OTS 44          1700  100 3.5  0.5 2300
KPNO Tau 4      1700  100 3.5  0.5 2300
2M0141          1800  200 4.5  0.5 2200
TWA 5B          2500  100 4.0  0.5 2570
Cha 1109        1800  100 4.0  0.5 2300
Gl 417 B        1800  100 5.0  0.5 1660
AB Pic b        1800  200 4.5  0.5 2200
CT Cha b        2600  100 3.5  0.5 2700
DH Tau B        2400  100 3.5  0.5 2350
GSC8047 B       2200  100 4.0  0.5 2300
USco CTIO 108B  2300  100 4.0  0.5 2300
HR7329 B        2600  100 . . . 2570
2M0345          2400  100 5.0  0.5 2230
TWA22 A         3000  100 4.5  0.5 3125
TWA22 B         3000  100 4.5  0.5 3064"""

path = "DI_spectra/"
files = ["L0_R2000_ABPicb_J_0.025_norm_master_spectrum.fits", "M8.5_R2000_2M1207A_SINFONIspeclib_JHK.fits", "M6_R2000_TWA8B.fits", "M5_R1200_TWA11C.fits"]

for i, f in enumerate(files):
    ax = plt.subplot(gs[i,1])
    p = pyfits.open(path + f)
    d = p[0].data 
    plt.plot(d[:,0], d[:,1])
    plt.xlim(1,2)


plt.show()
