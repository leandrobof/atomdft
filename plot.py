import numpy as np
import matplotlib.pyplot as plt

a=np.loadtxt("radial.txt")
b=np.loadtxt("radial.txt",usecols=(0,2))
plt.plot(a[:,0],a[:,1],b[:,0],b[:,1])
plt.show()
