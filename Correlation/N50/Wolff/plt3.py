import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

i=int(sys.argv[1])

a=pd.read_csv('cor.csv', header=None, sep='\s+') 
a=a.values
a=a[i-1,:]

plt.show()
#print(a)

s=0
for i in range(len(a)):
    if(a[i]<=0):
        break
    s+=a[i]
    #print(i)
print(s)
