from numpy import*
from matplotlib.pyplot import*
from celluloid import Camera

def f(u, mass, q1, q2, l):
    f = zeros(9)
    f[0] = u[4]
    f[1] = u[5]
    f[2] = u[6]
    f[3] = u[7]
    f[4] = q1*q2*u[0]/mass/(u[0]**2+u[1]**2)**(3/2)+2*u[8]/mass*(u[0]-u[2])
    f[5] = q1*q2*u[1]/mass/(u[0]**2+u[1]**2)**(3/2)+2*u[8]/mass*(u[1]-u[3])
    f[6] = -q1*q2*u[2]/mass/(u[2]**2+u[3]**2)**(3/2)-2*u[8]/mass*(u[0]-u[2])
    f[7] = q1*q2*u[3]/mass/(u[2]**2+u[3]**2)**(3/2)-2*u[8]/mass*(u[1]-u[3])
    f[8] = (u[0]-u[2])**2+(u[1]-u[3])**2-l**2
    return f

def Jacobian(u, mass):
    j = np.zeros((9,9))
    j[0,4] = 1.
    j[1,5] = 1.
    j[2,6] = 1.
    j[3,7] = 1.
    j[4,0] = q1*q2/mass*(u[1]**2-2*u[0]**2)/(u[0]**2+u[1]**2)**(5/2)+2*u[8]/mass
    j[4,1] = -3*q1*q2/mass*u[1]*u[0]/(u[0]**2+u[1]**2)**(5/2)
    j[4,2] = -2*u[8]/mass
    j[4,8] = 2*mass*(u[0]-u[2])
    j[5,0] = -3*q1*q2/mass*u[1]*u[0]/(u[0]**2+u[1]**2)**(5/2)
    j[5,1] = q1*q2/mass*(u[0]**2-2*u[1]**2)/(u[0]**2+u[1]**2)**(5/2)+2*u[8]/mass
    j[5,3] = -2*u[8]/mass
    j[5,8] = 2*mass*(u[1]-u[3])
    j[6,0] = -2*u[8]/mass
    j[6,2] = -q1*q2/mass*(u[3]**2-2*u[2]**2)/(u[2]**2+u[3]**2)**(5/2)+2*u[8]/mass
    j[6,3] = 3*q1*q2/mass*u[2]*u[3]/(u[2]**2+u[3]**2)**(5/2)
    j[6,8] = -2*mass*(u[0]-u[2])
    j[7,1] = -2*u[8]/mass
    j[7,2] = -3*q1*q2/mass*u[2]*u[3]/(u[2]**2+u[3]**2)**(5/2)
    j[7,3] = q1*q2/mass*(u[2]**2-2*u[3]**2)/(u[2]**2+u[3]**2)**(5/2)+2*u[8]/mass
    j[7,8] = -2*mass*(u[1]-u[3])
    j[8,0] = 2*(u[0]-u[2])
    j[8,1] = 2*(u[1]-u[3])
    j[8,2] = -2*(u[0]-u[2])
    j[8,3] = -2*(u[1]-u[3])
    return j

M = 500
t_0 = 0.; T = 23
tau = (T-t_0)/M
t = linspace(t_0, T, M+1)
x1 = 6; x2 = x1
y1 = 0.1; y2 = -0.1
mass = 1
l = sqrt((x1-x2)**2+(y1-y2)**2)
q1 = 10
q2 = 0.8


u = zeros((M+1, 9))
u[0,0] = x1
u[0,1] = y1
u[0,2] = x2
u[0,3] = y2
u[0,4] = 0.
u[0,5] = 0.
u[0,6] = 0.
u[0,7] = 0.
u[0,8] = (-(x1-x2)/mass*q1*q2*(x1+x2)/(x1**2+y1**2)**(3/2)-(y1-y2)**2/mass*q1*q2/(x2**2+y2**2)**(3/2))/(4*(x1-x2)**2 + 4*(y1-y2)**2)

D = eye(9)
D[8,8] = 0
a1 = 1

fig, ax = subplots(figsize=(7,7))
camera = Camera(fig)

for m in range(M):
    w1 = np.linalg.solve((D-a1*tau*Jacobian(u[m], mass)), f(u[m], mass, q1,q2,l))
    u[m+1] = u[m]+tau*real(w1)
    
    if (m+1)%10 == 0 : 
        ax.plot([u[m+1,0],u[m+1,2]],[u[m+1,1],u[m+1,3]],color='black')
        circle2 = Circle((u[m+1,0], u[m+1,1]), 0.05, color='black')
        circle3 = Circle((u[m+1,2], u[m+1,3]), 0.05, color='black')
        circle4 = Circle((0, 0), 0.2, color='black')
        ax.add_patch(circle2)
        ax.add_patch(circle3)
        ax.add_patch(circle4)
        camera.snap()
        
ax.grid()
animation = camera.animate()  
animation.save('oscillator.gif', writer = 'pillow')
xlim((-1,x1+1)); ylim((-(y1+1),y1+1))
show()

#plot(t,u[:,8])
