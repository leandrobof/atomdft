import numpy as np
import matplotlib.pyplot as plt

pi=3.14159265359
a=np.loadtxt("../mgo_band.dat")
b=np.loadtxt("mgo_band")
r=range(0,len(a[:,0]))
for i in range(5,len(a[0,:])):
    plt.plot(r,a[:,i],color="b")
    plt.plot(b[:,0],b[:,i],color="r")

plt.show()
