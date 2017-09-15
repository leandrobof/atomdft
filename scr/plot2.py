import numpy as np
import matplotlib.pyplot as plt

pi=3.14159265359
a=np.loadtxt("H.txt")

plt.plot(a[:,0],a[:,1],label="ss")
plt.plot(a[:,0],a[:,2],label="sp")
plt.plot(a[:,0],a[:,3],label="pp_s")
plt.plot(a[:,0],a[:,4],label="pp_pi")
plt.legend()
plt.show()
