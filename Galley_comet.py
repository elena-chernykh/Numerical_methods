from numpy import*
from matplotlib.pyplot import*

def f(u, G, Mass, mass):
    f = zeros(4)
    f[0] = u[2]
    f[1] = u[3]
    f[2] = -G*Mass/(u[0]**2+u[1]**2)*u[0]/sqrt(u[1]**2+u[0]**2) #cos(atan)
    f[3] = -G*Mass/(u[0]**2+u[1]**2)*u[1]/sqrt(u[1]**2+u[0]**2) #sin(atan)
    return f

def RK_1(u, tau, G, Mass, mass):
    b_1 = 1
    for m in range(M):
        w_1 = f(u[m], G, Mass, mass)
        u[m + 1] = u[m] + tau*b_1*w_1
        #u[m+1] = u[m]+tau*f(u[m],G,Mass, mass)
    return u

def RK_2(u, tau, G, Mass, mass):
    a_21 = 2/3
    b_1 = 1/4; b_2 = 3/4
    for m in range(M):
        w_1 = f(u[m], G, Mass, mass)
        w_2 = f(u[m]+tau*a_21*w_1, G, Mass, mass)
        u[m + 1] = u[m] + tau*(b_1*w_1+b_2*w_2)
    return u

def RK_3(u, tau, G, Mass, mass):
    a_21 = 1/2; a_31 = 0; a_32 = 3/4
    b_1 = 2/9; b_2 = 3/9; b_3 = 4/9
    for m in range(M):
        w_1 = f(u[m], G, Mass, mass)
        w_2 = f(u[m]+tau*a_21*w_1, G, Mass, mass)
        w_3 = f(u[m]+tau*(a_31*w_1+a_32*w_2), G, Mass, mass)
        u[m + 1] = u[m] + tau*(b_1*w_1+b_2*w_2+b_3*w_3)
    return u

def RK_4(u, tau, G, Mass, mass):
    a_21 = 1/2; a_31 = 0; a_32 = 1/2; a_41 = 0; a_42 = 0; a_43 = 1
    b_1 = 1/6; b_2 = 1/3; b_3 = 1/3; b_4 = 1/6
    for m in range(M):
        w_1 = f(u[m], G, Mass, mass)
        w_2 = f(u[m]+tau*a_21*w_1, G, Mass, mass)
        w_3 = f(u[m]+tau*(a_31*w_1+a_32*w_2), G, Mass, mass)
        w_4 = f(u[m]+tau*(a_41*w_1+a_42*w_2+a_43*w_3), G, Mass, mass)
        u[m + 1] = u[m] + tau*(b_1*w_1+b_2*w_2+b_3*w_3+w_4*b_4)
    return u

t_0 = 0.; T = 75.3*31557600
x_0 = 5248240000.; y_0 = 0.
v_0_x = 0; v_0_y = 0.91158 # в афелии
Mass = 1.99*10**30 #масса Солнца
mass = 0
G = 6.67*10**(-20) #гравит. постоянная

M = 10000
tau = (T-t_0)/M
t = linspace(t_0, T, M+1)

u = zeros((M+1, 4))
u[0,0] = x_0
u[0,1] = y_0
u[0,2] = v_0_x
u[0,3] = v_0_y

u = RK_4(u, tau, G, Mass, mass)

figure()
plot(u[:,0], u[:,1], '-ro', markersize=5)

xlim((-548240000,5500000000)); ylim((-1624120000,1624120000))
gca().set_aspect('equal', adjustable='box')

