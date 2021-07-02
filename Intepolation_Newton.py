import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from sympy.abc import x

N = 5
#x = [4.,5.,3., 1., 2.]
#u = [54.6, 148.4, 20.1, 2.7, 7.4]
x = [1.,2.,3., 4., 5., 6., 7., 8.]
u = [2.7, 7.4, 20.1, 54.6, 148.4, 403.4, 1096.6, 2981.]

U = np.zeros((N,N))
for j in range(N):
    U[j][0] = u[j]

for j in range(1, N):
    for i in range(N-j):
        U[i][j] = j*(U[i][j-1]-U[i+1][j-1])/(x[i]-x[i+j])

coord = 7.8
counter = 0
P = 0.
for i in range(N):
    mult = 1.
    for j in range(i):
        mult *= (coord - x[j])/float(j+1)
    P += U[0][i]*mult
print(P)
plt.plot(x, u, "or")
plt.plot(coord, P, "og")
plt.show()