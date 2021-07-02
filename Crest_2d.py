from numpy import * 
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import *
from matplotlib import cm

def plots(X,Y,w,u,j):
    fig = figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, w[:,:,j], cmap= 'inferno',rstride=1,cstride=1,linewidth=0.1)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('w')
    ax.set_title('Численное решение')
    show()
    fig = figure()
    ax2 = fig.add_subplot(111, projection='3d')
    ax2.plot_surface(X, Y, u[:,:,j], cmap= 'inferno',rstride=1,cstride=1,linewidth=0.1)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('u')
    ax2.set_title('Аналитическое решение')
    show()
    fig = figure()
    ax3 = fig.add_subplot(111, projection='3d')
    ax3.plot_surface(X, Y, w[:,:,j]-u[:,:,j], cmap= 'inferno',rstride=1,cstride=1,linewidth=0.1)
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_zlabel('err')
    ax3.set_title('Ошибка')
    show()

a=1 
b=1 
N=101
M=101
T=1

h_x = a/(N-1/2)
h_y = b/(M-1/2)
J=int(100+round(2*T*sqrt(1/h_x**2+1/h_y**2)))
tau=T/(J-1)

x = np.linspace(-1/2*h_x,a,N+1)
y = np.linspace(-1/2*h_y,b,M+1)
t = np.linspace(0,T,J)
u = np.zeros((N+1,M+1,J))
w = np.zeros((N+1,M+1,J))
X = np.zeros((N+1,M+1))
Y = np.zeros((N+1,M+1))
w1 = np.zeros((N+1,M+1))

# аналитическое решение ------------------------------------------------------
for n in range(N+1):
    for m in range(M+1):
        for j in range(J):
            u[n,m,j] = (t[j]**2) * (x[n]+y[m]) + 5* cos(pi*sqrt(10)*t[j]/2)* \
            cos(pi*x[n]/2)* cos(3*pi*y[m] /2)
for n in range(N+1):
    for m in range(M+1):
        X[n,m]=x[n]
        Y[n,m]=y[m]

# численное решение ----------------------------------------------------------
for n in range(N+1):
    for m in range(M+1):
        w[n,m,0] = 5*cos(pi*x[n]/2)*cos(3*pi*y[m]/2)
        w[n,m,1]=w[n,m,0] + tau**2/2*(-25*(pi**2)/4*np.cos(pi*x[n]/2)*cos(3*pi*y[m]/2)) + tau**2*(x[n]+y[m])
for n in range(N+1):
    w[n,M,1] = t[1]**2 * (x[n]+1)
    
for m in range(M+1):
    w[N,m,1] = t[1]**2 * (y[m]+1)
    
for j in range(1,J-1):
    for n in range(1,N):
        for m in range(1,M):
            w[n,m,j+1] = 2*w[n,m,j]  - w[n,m,j-1] + \
                (tau**2)*((w[n+1,m,j] - 2*w[n,m,j] + w[n-1,m,j])/(h_x**2) + \
                          (w[n,m+1,j] - 2*w[n,m,j] + w[n,m-1,j])/(h_y**2) +2*(x[n]+y[m]))
                
    for m in range(M+1):
        w[0,m,j+1] = w[1,m,j+1] - t[j+1]**2*h_x
        w[N,m,j+1] = t[j+1]**2*(1+y[m])

    for n in range(N+1):
        w[n,0,j+1] = w[n,1,j+1] - t[j+1]**2*h_y
        w[n,M,j+1] = t[j+1]**2*(1+x[n])

plots(X,Y,w,u,245)
'''err = zeros((N+1,M+1,J))
for j in range(J):
    for n in range(N+1):
        for m in range(M+1):
            err[n,m,j] = u[n,m,j]-w[n,m,j]
            
def data(i):
    ax.clear()
    surf = ax.plot_surface(y, x, err[:,:,i], cmap= 'inferno', linewidth=0, antialiased=False)
    ax.set_zlim([-0.06,0.06]) # set zlim to be always the same for every frame

fig = figure()
ax = fig.gca(projection='3d')
x,y = np.meshgrid(y, x)

ani = FuncAnimation(fig, data, J, interval=5, repeat=True )
title('Numerical solution')
xlabel("x")
ylabel("y")
show()'''
