# экономично-факторизованная схема------------------------------------------ 
from numpy import * 
from matplotlib.pyplot import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import *

L = 1
H = 1
T = 1
N = 51
M = 51
J = 121

x = zeros(N)
t = zeros(J)
y = zeros(M)

hx = L/(N-1-1/2)
for n in range(N):
    x[n] = -hx/2+n*hx
hy = H/(M-1-1/2)
for m in range(M):
    y[m] = -hy/2+m*hy

tau = T/(J-1)
t = linspace(0,T,J)
sigma = 0.25
#аналит.решение -------------------------------------------------------------
u = zeros((N,M,J))
max_u = zeros(J)
temp = zeros(N)

for j in range(J):
    for n in range(N):
        for m in range(M):
            u[n,m,j] = (t[j]**2)*(x[n]+y[m]) + 5*cos(pi*sqrt(10)/2*t[j])* \
                cos(pi*x[n]/2)*cos(3*pi*y[m]/2)
    for n in range(N):
        temp[n] = max(u[n,:,j])
    max_u[j] = max(temp)
    

w = zeros((N,M,J))
err = zeros((N,M,J))

for n in range(N):
    for m in range(M):
        w[n,m,0] = 5*cos(pi*x[n]/2)*cos(3*pi*y[m]/2)
        w[n,m,1]= w[n,m,0] - tau**2/2*(-25*(pi**2)/4*cos(pi*x[n]/2)*cos(3*pi*y[m]/2)) + tau**2*(x[n]+y[m])
for n in range(N):
    w[n,M-1,0] = t[0]**2 * (x[n]+1)  
for m in range(M):
    w[N-1,m,0] = t[0]**2 * (y[m]+1)

w1 = zeros((N,M))

Ax = tau**2*sigma/(hx**2)
Cx = 1+2*Ax
Ay = tau**2*sigma/(hy**2)
Cy = 2*Ay+1

alpha_x = zeros(N-1)
beta_x = zeros(N-1)
alpha_y = zeros(M-1)
beta_y = zeros(M-1)

for j in range(J-1):
    
    for m in range(1,M-1):
        alpha_x[0] = 1
        beta_x[0] = -2*hx
        for n in range(1,N-1):
            F = (w[n-1,m,j]-2*w[n,m,j]+w[n+1,m,j])/(hx**2) + \
                (w[n,m+1,j]-2*w[n,m,j]+w[n,m-1,j])/(hy**2) + 2*(y[m]+x[n])
            alpha_x[n] = Ax/(Cx-Ax*alpha_x[n-1])
            beta_x[n] = (F+Ax*beta_x[n-1])/(Cx-Ax*alpha_x[n-1])
        w1[N-1,m] = 2*(1+y[m])    
        for n in range(N-2,-1,-1):
            w1[n,m] = alpha_x[n]*w1[n+1,m]+beta_x[n]  
            
    for n in range(1,N-1):
        alpha_y[0] = 1
        beta_y[0] = -(t[j]**2)*hy
        for m in range(1,M-1):
            #F = tau**2*w1[n,m] + sigma*tau**2/(hy**2)*(w[n,m-1,j-1]-2*w[n,m-1,j]) - \
            #    (1+2*sigma*tau**2/(hy**2))*(w[n,m,j-1]-2*w[n,m,j])+ \
            #        sigma*tau**2/hy**2*(w[n,m+1,j-1]-2*w[n,m+1,j])
            F= tau**2 * w1[n,m] + Ay*(w[n,m-1,j-1]-2*w[n,m-1,j]) - \
                Cy*(w[n,m,j-1]-2*w[n,m,j]) + Ay*(w[n,m+1,j-1]- 2*w[n,m+1,j])
            alpha_y[m] = Ay/(Cy-Ay*alpha_y[m-1])
            beta_y[m] = (F+Ay*beta_y[m-1])/(Cy-Ay*alpha_y[m-1])
        w[n,M-1,j+1] = (t[j+1]**2)*(x[n]+1)
        err[n,M-1,j+1] = abs(u[n,M-1,j+1]-w[n,M-1,j+1])
        
        for m in range(M-2,-1,-1):
             w[n,m,j+1] = alpha_y[m]*w[n,m+1,j+1]+beta_y[m]
             err[n,m,j+1] = abs(u[n,m,j+1]-w[n,m,j+1])
    
    for m in range(M):
        w[0,m,j+1] = w[1,m,j+1]-hx*2*t[j+1]**2
        w[N-1,m,j+1] = t[j+1]**2*(y[m]+1)
        err[0,m,j+1] = abs(w[0,m,j+1]-u[0,m,j+1])
        err[N-1,m,j+1] = abs(w[N-1,m,j+1]-u[N-1,m,j+1])


fig = figure()
ax = fig.gca(projection='3d')
x,y = np.meshgrid(y, x)
surf = ax.plot_surface(y, x, u[:,:,0], cmap= 'inferno', linewidth=0, antialiased=False)
title('Numerical solution')
xlabel("x")
ylabel("y")
show()

'''def data(i):
    ax.clear()
    surf = ax.plot_surface(y, x, u[:,:,i], cmap= 'inferno', linewidth=0, antialiased=False)
    ax.set_zlim([-2,2]) # set zlim to be always the same for every frame

fig = figure()
ax = fig.gca(projection='3d')
x,y = np.meshgrid(y, x)

ani = FuncAnimation(fig, data, J, interval=10, repeat=True )
title('Analitical solution')
xlabel("x")
ylabel("y")
show()''' 


