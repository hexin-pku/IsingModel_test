import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

f1=str(sys.argv[1]) 
b0=float(sys.argv[2])

a=pd.read_csv(f1, sep='\s+')
a=a.values

m=a[:,0]

b=np.arange(len(m))*b0


plt.plot(b,m,'r--')
plt.ylabel('M [a.u.]')
plt.xlabel(r'$h$ [a.u.]')
plt.savefig('poltm.png')

plt.clf()
