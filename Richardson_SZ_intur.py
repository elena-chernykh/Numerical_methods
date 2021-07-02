import numpy as np
import sympy as sp

def make_Matrix(a, b, N, x, y):
    h = (b-a)/N
    A = np.zeros((N,N))
    for i in range(N):
        for j in range(1,N-1):
            A[i, j] = h*np.sin(x[i]+y[j])
    for i in range(N):
        A[i,0] = h*np.sin(x[i]+y[0])/2
        A[i,N-1] = h*np.sin(x[i]+y[N-1])/2
    return A

def Lambda(a, b, N, x, y):
    X_0 = np.ones(N)
    X_1 = np.ones(N)
    A = make_Matrix(a, b, N, x, y)
    for i in range(30):
        X_1 = np.dot(np.linalg.inv(A), X_0)
        l = np.dot(X_1, X_1)/np.dot(X_0, X_1)
        X_0 = X_1
    return l
    #while (abs(temp-l)/l > e):
    #    temp = l
    #    X_1 = np.dot(np.linalg.inv(A), X_0)
    #    l = np.dot(X_1, X_1)/np.dot(X_0, X_1)
    #    X_0 = X_1 
    #return 1/l

def PrintTriangular(A,i,s):
    for m in range(s):
        for n in range(m + 1 - i): 
            
            print('{0:7.4f}'.format(A[m,n]),end=' ')
        print()
    print()
 
a = 0.
b = np.pi
N = 2
r = 2
p = [1]
q = 1
S = 10
s = 0
eps = 0.03
check = True

for i in range(1,S):
    p.append(p[i-1]+q)

U = np.zeros((S, S))
R = np.zeros((S,S))
p_eff = np.zeros((S,S))

x = np.linspace(a,b, N)
y = np.linspace(a,b, N)
U[0,0] = Lambda(a, b, N, x, y)
i = 0
#for i in range(1, S):
while(check):
    i += 1
    x = np.linspace(a, b, N*r**i)
    y = np.linspace(a, b, N*r**i)
    U[i,0] = Lambda(a, b, N*r**i, x, y)
    for j in range(i):
        R[i, j] = (U[i, j] - U[i-1, j])/(r**p[j]-1)
        if (i >= 2 and R[i-1, j] != 0):
            p_eff[i, j] = np.log(abs(R[i-1, j]/R[i, j]))/np.log(r)
        U[i, j + 1] = R[i, j]+ U[i, j]
    s += 1
    if (s == S-1 or abs(R[i, j]/U[i, j+1]) < eps):
        check = False
        break

#print(s)
PrintTriangular(U,0,s)
PrintTriangular(R,1,s)
PrintTriangular(p_eff,2,s)
