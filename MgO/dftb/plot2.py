import numpy as np
import matplotlib.pyplot as plt
import sys

pi=3.14159265359
a=np.loadtxt("../mgo_band.dat")
b=np.loadtxt("mgo_band_tot.dat")
r=range(0,len(b[:,0]))
Efermi=sys.argv[1]
efermi=float(Efermi)*np.ones(np.shape(a))
b=np.sort(b)

for i in range(4,len(a[0,:])):
    plt.plot(r,a[:,i],color="r")

for j in range(0,len(b[0,:])):
    plt.plot(r,(b[:,j]-efermi[:,1]),color="b")

plt.show()
