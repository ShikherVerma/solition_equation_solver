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


# Write code for solving dy_dt using RK4 method

#def func(t, y):
    #"""
    #"""
    #return np.sin(t)

#def y_anal(t):
    #"""
    #"""
    #return -np.cos(t)

#def f(t, uHat, k):
    #"""
    #"""
    #u = np.fft.ifft(uHat)
    #nonlinear_term = ((2 * np.pi * 1j * k) / 2) * uHat
    #dissipation_term = (2 * np.pi * 1j)**3 * k * uHat * 1.
    #rhs = dissipation_term - nonlinear_term
    #return rhs

dt = 1e-7
t = np.arange(0, 0.08, dt)
N_t = t.shape[0]

dx = 0.006
x = np.arange(0, 1, dx)
N_x = x.shape[0]
u = np.zeros([N_t, N_x])
# uHat = np.zeros([N_t, N_x])

x_c = 0.5
sigma = 0.1
u[0] = np.exp(-(x - x_c)**2 / sigma**2)
# u[0] = np.sin(np.pi * x)

uHat = np.fft.fft(u[0])
k    = np.fft.fftfreq(N_x, d=(1 / N_x))

energy = [np.trapz(u[0]**2, x, dx)]
# k1 = f(t[0], uHat, k)
# plt.plot(x, k1)
for t_idx in trange(N_t - 1, leave = False):
    uHat  = np.fft.fft(u[t_idx])
    uHat2 = np.fft.fft(u[t_idx]**2)
    
    u[t_idx + 1] = np.fft.ifft((-1j * 1 * np.pi * k * dt * uHat2 + uHat) / (1 - dt * k * 0.1 * (2 * np.pi * 1j)**3)).real
    #plt.clf()
    #if (t_idx % 10 == 0):
        #plt.plot(x, u[t_idx])
        #plt.ylim([0, 1.0101])
        #plt.xlim([0, 1.])
        #plt.xlabel("$x$")
        #plt.ylabel("$u$")
        ##plt.title('$\mathrm{Time}\; \mathrm{evolution}\; \mathrm{of}\; u = e^{(x - 0.5)^2 / 0.025^2}$')
        #plt.title(r'Time = ' + str(t[t_idx]))
        #plt.savefig('rho' + '%05d'%t_idx + '.png')
        ##plt.show()
        #plt.clf()
    #print("plot {}".format(t_idx))
    energy.append(np.trapz(u[t_idx + 1]**2, x, dx))
    pass

x1,x2,y1,y2 = plt.axis()
#plt.axis((0.0,0.08,0.12,0.13))
plt.plot(t, np.abs(np.array(energy)))
#plt.savefig('Energy1e-7.png')
plt.show()
print("done")
# plt.plot(t, np.abs(y - y_anal(t)))
# plt.plot(t, y_anal(t))