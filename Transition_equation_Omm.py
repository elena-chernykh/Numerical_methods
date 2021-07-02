from numpy import *
from matplotlib.pyplot import *
from matplotlib.animation import *
from mpl_toolkits.mplot3d import Axes3D

def u_init(x):
    u_init = x**2+1
    return u_init

def u_left(t):
    u_left = exp(-t)
    return u_left

def f(y,t,h,u_left):
    f = zeros(N)
    f[0] = -2*y[0]/(1+(1+y[0]**2)**2)*(y[0] - u_left(t))/h
    for n in range(1,N):
        f[n] = -2*y[n]/(1+(1+y[n]**2)**2)*(y[n] - y[n - 1])/h
    return f

def f_y(y,t,h,u_left):
    f_y = zeros((N,N))
    f_y[0,0] = 1/h*(u_left(t)*(-6*y[0]**4-4*y[0]**2+4)-4*y[0]*(y[0]**4-2))/(y[0]**4+2*y[0]**2+2)**2
    for n in range(1,N):
        f_y[n,n] = 1/h*(y[n - 1]*(-6*y[n]**4-4*y[n]**2+4)-4*y[n]*(y[n]**4-2))/(y[n]**4+2*y[n]**2+2)**2
        f_y[n,n - 1] = -1/h*2*y[n]/(1+(1+y[n]**2)**2)
    return f_y

a = 0.; b = -1.
t_0 = 0.; T = 1.

N = 200
M = 300

a_11 = (1 + 1j)/2

h = abs(b - a)/N
x = linspace(a,b,N + 1)

tau = (T - t_0)/M
t = linspace(t_0,T,M + 1)

u = zeros((M + 1,N + 1))
y = zeros((M + 1,N))

for n in range(N + 1):
    u[0,n] = u_init(x[n])
    
y[0] = u[0,1:N+1]

for m in range(M):
    w_1 = linalg.solve(eye(N) - a_11*(t[m + 1] - t[m])*f_y(y[m],t[m],h,u_left),f(y[m],(t[m + 1] + t[m])/2,h,u_left))
    y[m + 1] = y[m] + (t[m + 1] - t[m])*w_1.real
    u[m + 1,0] = u_left(t[m])
    u[m + 1,1:N+1] = y[m + 1]
    
'''fig = figure() #построение графика
ax = fig.gca(projection = '3d')
x, t = meshgrid(x, t)
surf = ax.plot_surface(t,x,u, cmap = 'inferno')

xlabel("t")
ylabel("x")'''

fig = figure()
ax = axes(xlim=(a, b), ylim=(0, 2.0))
line, = ax.plot([], [], lw=3)

def init():
    line.set_data([], [])
    return line,
def animate(i):
    line.set_data(x,u[i,:])
    return line,

anim = FuncAnimation(fig, animate, init_func=init,frames = M + 1, interval = 50, blit=True)


