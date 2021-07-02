from numpy import*
from matplotlib.pyplot import*

def f(u, G, Mass):
    f = zeros(4)
    g_1 = -G*Mass*u[0]/(sqrt((u[0])**2+(u[1])**2))**3
    g_2 = -G*Mass*u[1]/(sqrt((u[0])**2+(u[1])**2))**3
    f[0] = u[2]/sqrt(1+u[2]**2+u[3]**2+g_1**2+g_2**2)
    f[1] = u[3]/sqrt(1+u[2]**2+u[3]**2+g_1**2+g_2**2)
    f[2] = g_1/sqrt(1+u[2]**2+u[3]**2+g_1**2+g_2**2)
    f[3] = g_2/sqrt(1+u[2]**2+u[3]**2+g_1**2+g_2**2)
    return f

def RK_1(u, dl, G, Mass, M):
    b_1 = 1
    for m in range(M):
        w_1 = f(u[m], G, Mass)
        u[m + 1] = u[m] + dl*b_1*w_1
    return u

def RK_2(u, dl, G, Mass, M):
    a_21 = 2/3
    b_1 = 1/4; b_2 = 3/4
    for m in range(M):
        w_1 = f(u[m], G, Mass)
        w_2 = f(u[m]+dl*a_21*w_1, G, Mass)
        u[m + 1] = u[m] + dl*(b_1*w_1+b_2*w_2)
    return u

def RK_3(u, dl, G, Mass, M):
    a_21 = 1/2; a_31 = 0; a_32 = 3/4
    b_1 = 2/9; b_2 = 3/9; b_3 = 4/9
    for m in range(M):
        w_1 = f(u[m], G, Mass)
        w_2 = f(u[m]+dl*a_21*w_1, G, Mass)
        w_3 = f(u[m]+dl*(a_31*w_1+a_32*w_2), G, Mass)
        u[m + 1] = u[m] + dl*(b_1*w_1+b_2*w_2+b_3*w_3)
    return u

def RK_4(u, dl, G, Mass, M):
    a_21 = 1/2; a_31 = 0; a_32 = 1/2; a_41 = 0; a_42 = 0; a_43 = 1
    b_1 = 1/6; b_2 = 1/3; b_3 = 1/3; b_4 = 1/6
    for m in range(M):
        w_1 = f(u[m], G, Mass)
        w_2 = f(u[m]+dl*a_21*w_1, G, Mass)
        w_3 = f(u[m]+dl*(a_31*w_1+a_32*w_2), G, Mass)
        w_4 = f(u[m]+dl*(a_41*w_1+a_42*w_2+a_43*w_3), G, Mass)
        u[m + 1] = u[m] + dl*(b_1*w_1+b_2*w_2+b_3*w_3+w_4*b_4)
    return u

def norm(r_err, dl, M): #норма
    norma = 0.
    for i in range(M):
        for j in range(4): # цикл по всем элементам
            norma += r_err[i,j]*r_err[i,j]*dl
    return sqrt(norma)

t_0 = 0.; T = 75.3*31557600
x_0 = 5248240000.; y_0 = 0.
v_0_x = 0; v_0_y = 0.91158 # в афелии
Mass = 1.99*10**30 #масса Солнца
mass = 0
G = 6.67*10**(-20) #гравит. постоянная

ex = 0.9671429
a = 2667950000
b = a*sqrt(1-ex**2) 
L = 4*(3.1415*a*b+(a-b)**2)/(a+b)*1.03
M = 300
dl_0 = L/M

u = zeros((M+1, 4))
u[0,0] = x_0
u[0,1] = y_0
u[0,2] = v_0_x
u[0,3] = v_0_y

r = 2
S = 8
p = 3

R = np.zeros(S)
p_eff = np.zeros(S)
percents = np.zeros(S)

u_1 = RK_4(u, dl_0, G, Mass, M)

for i in range(1,S):
    
    R_err = zeros((M+1, 4)) # для хранения погрешности
    dl = L/M/r
    u_2 = zeros((M*r+1, 4))
    u_2[0,0] = x_0; u_2[0,1] = y_0; u_2[0,2] = v_0_x; u_2[0,3] = v_0_y
    u_2 = RK_3(u_2, dl, G, Mass, M*r)
    
    for j in range(M):
        R_err[j] = (u_1[j]-u_2[j*r])/(r**p-1)
    R[i] = norm(R_err, dl_0, M)
    percents[i] = R[i]/norm(u_2, dl, M*r)*100
    if (i >= 2):
       p_eff[i] = log(R[i-1]/R[i])/log(r)
    u_1 = u_2
    M = M*r
    dl_0 = dl

for i in range(2,S):
    print("%3d %.5f %.5f" % (i+1, p_eff[i], percents[i]))
    
u = u_2
figure()
plot(u[:,0], u[:,1], '-ro', markersize=5)

xlim((-548240000,5500000000)); ylim((-1624120000,1624120000))
gca().set_aspect('equal', adjustable='box')


