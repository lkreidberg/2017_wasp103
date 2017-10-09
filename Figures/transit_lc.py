import batman
import numpy as np


def lc(t, rp): 
	params = batman.TransitParams()

	params.t0 = 2457080.64041702
	params.per =  0.925545613
	params.inc = 87.3
	params.ecc = 0. 
	params.w = 90. 
	params.rp = rp
	params.a = 3.0
	params.limb_dark = "quadratic"
	params.u = [0.1, 0.3]

	m = batman.TransitModel(params, t)
	lc = m.light_curve(params)
	
	return lc 
