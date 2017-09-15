import numpy as np
import matplotlib.pyplot as plt

pi=3.14159265359
b=np.loadtxt("radial2.txt")
for i in range(1,len(b[0,:])):
    plt.plot(b[:,0],b[:,i])
plt.show()
