import numpy as np
import sys
import math

archivo=open(sys.argv[1],"r")
datos=archivo.readline()
datos=datos.split()
nbd=datos[2]
nbd=nbd.split(',')
nbd=nbd[0]
nks=datos[4]
band=np.zeros((int(nks),int(nbd)))
Efermi=float(sys.argv[2])

for i in range(int(nks)):
   archivo.readline()
   for l in range(int(math.ceil(int(nbd)/10.))):
       k=archivo.readline()
       k=k.split()
       for j in range(len(k)):
            print (float(k[j])-Efermi) ,
       
   print 




archivo.close()
