import spiderman
import numpy as np

def proj_area(phi, inc):
        R0 = 1.07e8             #planet radius in m
        p = 0.00117             # planet-to-star mass ratio
        r = 2.9695e9            #orbital separation in m
	kn = 0.653		#from table B5 of Leconte et al. 2011
        n = 1.005
        qn = kn*(1.- n/5.)      #n 

        alpha1 = 2.5*qn*R0**3/(p*r**3)
        alpha2 = -1.25*qn*R0**3/(p*r**3)
        alpha3 = alpha2
        
	a1 = R0*(1+alpha1)
        a2 = R0*(1+alpha2)
        a3 = R0*(1+alpha3)

        return  np.sqrt(a3**2*np.sin(inc)**2*(a1**2*np.sin(phi)**2+a2**2*np.cos(phi)**2)+ a1**2*a2**2*np.cos(inc)**2)/(a2*a3)


def lc(t, rp, T_s, l1, l2, xi, T_n, delta_T, dilution, eclipse):
	web_p = spiderman.ModelParams(brightness_model =  'zhang') 

	web_p.n_layers = 5
	#web_p.t0 = 0. 
	web_p.t0 = 2457080.64041702
	web_p.per =  0.925545613
	web_p.a_abs = 1.3260e-02
	web_p.inc = 87.3
	web_p.ecc = 0. 
	web_p.w = 90. 
	web_p.rp = 0.115
	web_p.a = 3.0
	web_p.p_u1 = 0. 
	web_p.p_u2 = 0. 
	web_p.T_s = T_s 
	web_p.l1 = l1 
	web_p.l2 = l2 
	web_p.xi = xi 
	web_p.T_n = T_n 
	web_p.delta_T = delta_T
	web_p.delta_T_cloud = 0.
	web_p.thermal = True
        web_p.dilution = dilution
	web_p.eclipse = eclipse

	phs = (t - web_p.t0)/web_p.per
	phs -= np.round(phs)
	rprs2 = proj_area(phs*2.*np.pi, web_p.inc*np.pi/180.)
	lc = spiderman.web.lightcurve(t, web_p)
	
	return (np.array(lc) - 1.0)*rprs2 + 1.
