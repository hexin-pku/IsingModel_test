import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

N = int(sys.argv[1])
H = np.zeros(N)
C = np.zeros(N)
G = np.zeros(N)
X = np.zeros(N)
A = np.zeros((N,4))

for i in range(1,N+1):
    a = pd.read_csv('beta_i%d'%i+'.dat', header=None, sep='\s+')
    a = a.values
    H[i-1] = a[:,1].mean()
    C[i-1] = a[:,1].std()
    G[i-1] = a[:,2].mean()
    X[i-1] = a[:,2].std()

A[:,0] = H
A[:,1] = C
A[:,2] = G
A[:,3] = X

a = pd.DataFrame(A)
a.to_csv('HCGX.dat')

plt.plot(np.linspace(1,N,N), H, 'r--')
plt.savefig('H.png')
plt.clf()

plt.plot(np.linspace(1,N,N), C, 'r--')
plt.savefig('C.png')
plt.clf()

plt.plot(np.linspace(1,N,N), G, 'r--')
plt.savefig('G.png')
plt.clf()

plt.plot(np.linspace(1,N,N), X, 'r--')
plt.savefig('X.png')

