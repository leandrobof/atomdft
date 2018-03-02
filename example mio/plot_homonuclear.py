import numpy as np
import matplotlib.pyplot as plt


#a=np.loadtxt("C-C.skf")
b=np.loadtxt("H-H.skf")
#c=np.loadtxt("C-C2.skf",skiprows=22)
d=np.loadtxt("H-H2.skf",skiprows=22)

r2=np.linspace(0.4,12,581)
r1=np.linspace(0.4,10.4,500)

for i in range(0,20):
    #plt.plot(r2,c[:,i],color="r")     
    plt.plot(r2,d[:,i],color="r")


for j in range(0,len(b[1,:])):
     plt.plot(r1,b[:,j],color="b")




#plt.plot(r1,a[:,0],color="b")
#plt.plot(r1,a[:,1],color="b")     
#plt.plot(r1,a[:,2],color="b")
#plt.plot(r1,a[:,3],color="b")


plt.show()
