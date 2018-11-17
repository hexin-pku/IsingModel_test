import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

i=int(sys.argv[1])

a=pd.read_csv('cor.csv', header=None, sep='\s+') 
a=a.values
a=a[i,:]

plt.plot(np.arange(len(a)), a, 'r--')
plt.show()

print(a.sum())
