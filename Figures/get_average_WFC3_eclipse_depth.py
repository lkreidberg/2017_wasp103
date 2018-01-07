import numpy as np

def weighted_mean(data, err):            #calculates the weighted mean for data points data with std devs. err
	ind = err != 0.0
        weights = 1.0/err[ind]**2
        mu = np.sum(data[ind]*weights)/np.sum(weights)
        var = 1.0/np.sum(weights)
        return [mu, np.sqrt(var)]                


d = np.genfromtxt("espec_dayside.txt")

print weighted_mean(d[0:10,1], d[0:10,2])
