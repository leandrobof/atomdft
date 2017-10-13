import numpy as np
import matplotlib.pyplot as plt

pi=3.14159265359
a=b=np.loadtxt("radial.txt")
b=np.loadtxt("radial2.txt")
for i in range(1,len(b[0,:])):
    noconf,=plt.plot(a[:,0],a[:,i],color="b")

for i in range(1,len(b[0,:])):
    conf,=plt.plot(b[:,0],b[:,i],color="r")

plt.legend([noconf,conf],["no conf","conf"])
plt.show()
