import matplotlib.pyplot as plt
import numpy as np
import os, glob

gcm_path = "./Vivien_models2"
gcms = glob.glob(os.path.join(gcm_path, "PC*PT.dat"))


#plot gcms
for g in gcms:
    print g,
    model = np.genfromtxt(g, delimiter = ',') 
    plt.plot(model[:,0], model[:,5]*1e3)
    print (np.max(model[:,5]) - np.min(model[:,1]))*1e3
    

plt.show()
