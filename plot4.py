import numpy as np
import matplotlib.pyplot as plt

pi=3.14159265359
a=np.loadtxt("H.txt")
b=np.loadtxt("H2.txt")

plt.plot(a[:,0],a[:,1],color="b")
plt.plot(a[:,0],a[:,2],color="b")
plt.plot(a[:,0],a[:,3],color="b")
plt.plot(a[:,0],a[:,4],color="b")

plt.plot(b[:,0],b[:,1],color="r")
plt.plot(b[:,0],b[:,2],color="r")
plt.plot(b[:,0],b[:,3],color="r")
plt.plot(b[:,0],b[:,4],color="r")


plt.legend()
plt.show()
