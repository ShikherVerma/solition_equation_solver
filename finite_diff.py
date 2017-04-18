import matplotlib.pyplot as plt
import numpy as np
import math
from scipy import integrate
from tqdm import trange
from scipy import special as sp
import matplotlib.lines as lines
import inspect

plt.rcParams['figure.figsize']  = 12, 7.5
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['font.family']     = 'serif'
plt.rcParams['font.weight']     = 'bold'
plt.rcParams['font.size']       = 20  
plt.rcParams['font.sans-serif'] = 'serif'
plt.rcParams['text.usetex']     = True
plt.rcParams['axes.linewidth']  = 1.5
plt.rcParams['axes.titlesize']  = 'medium'
plt.rcParams['axes.labelsize']  = 'medium'

plt.rcParams['xtick.major.size'] = 8
plt.rcParams['xtick.minor.size'] = 4
plt.rcParams['xtick.major.pad']  = 8
plt.rcParams['xtick.minor.pad']  = 8
plt.rcParams['xtick.color']      = 'k'
plt.rcParams['xtick.labelsize']  = 'medium'
plt.rcParams['xtick.direction']  = 'in'    

plt.rcParams['ytick.major.size'] = 8
plt.rcParams['ytick.minor.size'] = 4
plt.rcParams['ytick.major.pad']  = 8
plt.rcParams['ytick.minor.pad']  = 8
plt.rcParams['ytick.color']      = 'k'
plt.rcParams['ytick.labelsize']  = 'medium'
plt.rcParams['ytick.direction']  = 'in'
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True

def d1_backward(h, u0, u1, u2):
    return ( 3*u2
            -4*u1
            +  u0 )/(2*h)

def d2_backward(h, u0, u1, u2, u3, u4):
    return ( 3*d1_backward(h, u2, u3, u4)
            -4*d1_backward(h, u1, u2, u3)
            +  d1_backward(h, u0, u1, u2) )/(2*h)

def d3_backward(h, u0, u1, u2, u3, u4, u5, u6):
    return ( 3*d2_backward(h, u2, u3, u4, u5, u6)
            -4*d2_backward(h, u1, u2, u3, u4, u5)
            +  d2_backward(h, u0, u1, u2, u3, u4) )/(2*h)

def d1(h, u0, u2):
    return (u2-u0)/(2*h)

def d2(h, u0, u1, u2):
    return (u2 - 2*u1 + u0)/h**2

def d3(h, u0, u1, u2, u3, u4):
    return d1(h, d2(h, u0, u1, u2), d2(h, u2, u3, u4))

def d1_forward(h, u0, u1, u2):
    return (-3*u0
            +4*u1
            -  u2 )/(2*h)

def d2_forward(h, u0, u1, u2, u3, u4):
    return (-3*d1_forward(h, u0, u1, u2)
            +4*d1_forward(h, u1, u2, u3)
            -  d1_forward(h, u2, u3, u4) )/(2*h)

def d3_forward(h, u0, u1, u2, u3, u4, u5, u6):
    return (-3*d2_forward(h, u0, u1, u2, u3, u4)
            +4*d2_forward(h, u1, u2, u3, u4, u5)
            -  d2_forward(h, u2, u3, u4, u5, u6) )/(2*h)

dt = 1e-5
t = np.arange(0, 0.08, dt)
N_t = t.shape[0]

dx = 0.006
x = np.arange(0, 1, dx)
N_x = x.shape[0]
u = np.zeros([N_t, N_x])

x_c = 0.5
sigma = 0.1
u[0] = np.exp(-(x - x_c)**2 / sigma**2)
neu = 1
energy = [np.trapz(u[0]**2, x, dx)]

for t_idx in trange(N_t - 1, leave = False):
    for i in range(N_x):
        der1 = 0;
        der3 = 0;
        if i >= 2 and i<=N_x-3:
            # central method
            der1 = d1(dx, u[t_idx][i-1]**2, u[t_idx][i+1]**2)
            der3 = d3(dx, u[t_idx][i-2], u[t_idx][i-1], u[t_idx][i], u[t_idx][i+1], u[t_idx][i+2])
        elif i<2:
            # forward method
            der1 = d1_forward(dx, u[t_idx][i]**2, u[t_idx][i+1]**2, u[t_idx][i+2]**2)
            der3 = d3_forward(dx, u[t_idx][i], u[t_idx][i+1], u[t_idx][i+2], u[t_idx][i+3], u[t_idx][i+4], u[t_idx][i+5], u[t_idx][i+6])
        elif i>N_x-3:
            # backward method
            der1 = d1_backward(dx, u[t_idx][i-2]**2, u[t_idx][i-1]**2, u[t_idx][i]**2)
            der3 = d3_backward(dx, u[t_idx][i-6], u[t_idx][i-5], u[t_idx][i-4], u[t_idx][i-3], u[t_idx][i-2], u[t_idx][i-1], u[t_idx][i])

        u[t_idx+1][i] = u[t_idx][i] - 0.5 * dt * der1 + neu * dt * der3
    
    if (t_idx % 10 == 0):
        plt.plot(x, u[t_idx])
        plt.ylim([0, 1.0101])
        plt.xlim([0, 1.])
        plt.xlabel("$x$")
        plt.ylabel("$u$")
        plt.title(r'Time = ' + str(t[t_idx]))
        plt.savefig('_rho' + '%05d'%t_idx + '.png')
        plt.clf()
        pass
    
    energy.append(np.trapz(u[t_idx + 1]**2, x, dx))
    pass

#x1,x2,y1,y2 = plt.axis()
#plt.axis((0.0,0.08,0.12,0.13))
plt.plot(t, np.abs(np.array(energy)))
#print(energy)
#plt.savefig('Energy1e-7.png')
plt.show()
print("done")
