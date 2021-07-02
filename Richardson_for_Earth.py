from numpy import*
from matplotlib.pyplot import*

def f(u, G, Mass):
    f = zeros(4)
    f[0] = u[2]
    f[1] = u[3]
    f[2] = -G*Mass/(u[0]**2+u[1]**2)*u[0]/sqrt(u[1]**2+u[0]**2) #cos(atan)
    f[3] = -G*Mass/(u[0]**2+u[1]**2)*u[1]/sqrt(u[1]**2+u[0]**2) #sin(atan)
    return f

'''def RK_2(u, tau, G, Mass, M):
    a_21 = 2/3
    b_1 = 1/4; b_2 = 3/4
    for m in range(M):
        w_1 = f(u[m], G, Mass)
        w_2 = f(u[m]+tau*a_21*w_1, G, Mass)
        u[m + 1] = u[m] + tau*(b_1*w_1+b_2*w_2)
    return u'''
    
def RK_4(u, tau, G, Mass, M):
    a_21 = 1/2; a_31 = 0; a_32 = 1/2; a_41 = 0; a_42 = 0; a_43 = 1
    b_1 = 1/6; b_2 = 1/3; b_3 = 1/3; b_4 = 1/6
    for m in range(M):
        w_1 = f(u[m], G, Mass)
        w_2 = f(u[m]+tau*a_21*w_1, G, Mass)
        w_3 = f(u[m]+tau*(a_31*w_1+a_32*w_2), G, Mass)
        w_4 = f(u[m]+tau*(a_41*w_1+a_42*w_2+a_43*w_3), G, Mass)
        u[m + 1] = u[m] + tau*(b_1*w_1+b_2*w_2+b_3*w_3+w_4*b_4)
    return u

def norm(r_err, tau, M): #норма
    norma = 0.
    for i in range(M):
        for j in range(4): # цикл по всем элементам
            norma += r_err[i,j]*r_err[i,j]*tau
    return sqrt(norma)
        
        

t_0 = 0.; T = 31556925.
x_0 = 152000000.; y_0 = 0.
v_0_x = 0; v_0_y = 29.4 # в афелии
mass = 5.97*10**24 #масса Земли
Mass = 1.99*10**30 #масса Солнца
G = 6.67*10**(-20) #гравит. постоянная

M = 10
tau_0 = (T-t_0)/M
t = linspace(t_0, T, M+1)

r = 2
S = 8
p = 4


u = zeros((M+1, 4))
u[0,0] = x_0
u[0,1] = y_0
u[0,2] = v_0_x
u[0,3] = v_0_y

R = np.zeros(S)
p_eff = np.zeros(S)
percents = np.zeros(S)

u_1 = RK_4(u, tau_0, G, Mass, M)

for i in range(1,S):
    
    R_err = zeros((M+1, 4)) # для хранения погрешности
    tau = (T-t_0)/M/r
    u_2 = zeros((M*r+1, 4))
    u_2[0,0] = x_0; u_2[0,1] = y_0; u_2[0,2] = v_0_x; u_2[0,3] = v_0_y
    u_2 = RK_4(u_2, tau, G, Mass, M*r)
    
    for j in range(M):
        R_err[j] = (u_1[j]-u_2[j*r])/(r**p-1)
    R[i] = norm(R_err, tau_0, M)
    percents[i] = R[i]/norm(u_2, tau, M*r)*100
    if (i >= 2):
       p_eff[i] = log(R[i-1]/R[i])/log(r)
    u_1 = u_2
    M = M*r
    tau_0 = tau

for i in range(2,S):
    print("%3d %.5f %.5f" % (i+1, p_eff[i], percents[i]))
    
u = u_2   
figure()
plot(u[:,0], u[:,1], '-ro', markersize=5)

xlim((-162000000,162000000)); ylim((-162000000,162000000))
gca().set_aspect('equal', adjustable='box')

