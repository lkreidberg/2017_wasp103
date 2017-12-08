import matplotlib.pyplot as plt
import numpy as np
import os, glob

gcm_path = "./GCM_From_Vivien/GCMmodels"
gcms = glob.glob(os.path.join(gcm_path, "PC*PT.dat"))


#plot gcms
for g in gcms:
    model = np.genfromtxt(g, delimiter = ',') 
    plt.plot(model[:,0]/np.max(model[:,0])/2. + 0.5, model[:,5]*1e3)
    

plt.show()
