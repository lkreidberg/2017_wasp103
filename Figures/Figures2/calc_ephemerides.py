import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

def residuals(params, tt, tt_err):
	per, ephem = params
	t_num = np.round((tt - ephem)/per)
	model = per*t_num
        print "chisq", np.sum(((tt - ephem - model)/tt_err)**2)
        return (tt - ephem - model)/tt_err

#Measured transit times (time of inferior conjunction)
# 2456459.59957 +/- 0.00075	; Gillon 2014
# 2456836.296445 +/- 0.000055   ; Southworth 2015 (includes Gillon point)

#Ch2:
#Eclipse 1 midpoint: 	2457162.54885	+/- 0.00095
#Eclipse 2 midpoint: 	2457163.4775	+/- 0.0011
#Mid-transit time:	2457163.01496	+/- 0.00034
                    
#Ch1:
#Eclipse 1 midpoint: 	2457170.8814	+/- 0.0014	
#Eclipse 2 midpoint:	2457171.8095	+/- 0.0012	 
#Mid-transit time:	2457171.34370	+/- 0.00035 

#transit times and uncertainties (does not use Gillon point because Southworth et al. include it)
#tt = np.array([2456836.296445, 2457163.01496, 2457171.34370])
#tt_err = np.array([0.000055, 0.00034, 0.00035])


tt = np.array([])


per0 = 0.925545613
ephem0 = tt[0]

p0 = [per0, ephem0]
plsq, success  = leastsq(residuals, p0, args=(tt, tt_err))
per, ephem = plsq
#print "Not using best fit"
#per, ephem = per0, ephem0
print "difference in period", (per - per0)*24.*60.*60.

t_num = np.round((tt - ephem)/per)
tt_pred = per*t_num

plt.axhline(0, linestyle='dashed', color='0.5')
plt.errorbar(t_num, (tt - ephem - tt_pred)*24.*60., yerr = tt_err*24.*60., linestyle = 'none', marker = 'o')
#plt.xlabel("Cycle")
plt.ylabel("O - C (minutes)")
#plt.xlim((-10, 400))
plt.xlim((350, 365))


plt.show()


