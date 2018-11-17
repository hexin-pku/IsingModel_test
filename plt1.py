import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

f1=str(sys.argv[1]) 
b0=float(sys.argv[2])

a=pd.read_csv(f1, sep='\s+')
a=a.values
h=a[:,0]
c=a[:,1]
m=a[:,2]
x=a[:,3]
b=np.arange(len(h))*b0

plt.plot(b,h,'r--')
plt.ylabel('H [a.u.]')
plt.xlabel(r'$\beta$')
plt.savefig('polth.png')
plt.clf()

plt.plot(b,c,'r--')
plt.ylabel('C [a.u.]')
plt.xlabel(r'$\beta$')
plt.savefig('poltc.png')
plt.clf()

plt.plot(b,m,'r--')
plt.ylabel('M [a.u.]')
plt.xlabel(r'$\beta$')
plt.savefig('poltm.png')
plt.clf()

plt.plot(b,x,'r--')
plt.ylabel('X [a.u.]')
plt.xlabel(r'$\beta$')
plt.savefig('poltx.png')
plt.clf()
