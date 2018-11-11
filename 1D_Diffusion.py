# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 16:10:54 2016

@author: bennski

-----------------------------------------------------------------------------
 Backward time, centered space (BTCS) scheme to solve the 1D diffusion equation 
 for a molecules outflowing from a meteorite into a warm little pond.

 The equation solved is

       du     d  du
   Phi -- = D -- --
       dt     dx dx
       
 Phi:       porosity
 D:         diffusion coefficient
 x = 0:     center of meteorite
 x = xmax:  edge of meteorite
-----------------------------------------------------------------------------
"""

from matplotlib import pyplot as plt
import matplotlib.animation as animation
import numpy as np


# Initialize writer
Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

# Initialization function for animation
def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

# Set the diffusion of mass coefficient
D = (5.36e-11/0.25)*3600.*24. # diffusion/porosity [m^2 day^-1]
T = 273.15+65

xmax = 1e-2 # 1 cm radius meteorite
tmax = 10 # days
level = 8

nx = (2**level) + 1
nt = (2**(level+1)) + 1 # choosing nt to have twice as many grid points as nx

# Concentration array
u = np.matlib.zeros(shape=(nt, nx))

# Arrays for plotting
x = np.linspace(0,xmax,nx)
t = np.linspace(0,tmax,nt)

# Calculate delta_x, delta_t
delta_x = x[1] - x[0]
delta_t = t[2] - t[1]

# Initial condition
maxconc = 234
smooth = 1.
sum_init = 0 # sum of the initial concentration
for k in xrange(0,nx-1):
    u[0,k] = maxconc*(1-np.e**(-smooth*(nx-1-k))) # smoothed outer edge
    sum_init = sum_init + u[0,k]
    
# For plotting initial condition
r1 = np.zeros(shape=(nx))
for r in xrange(0,nx):
    r1[r] = u[0,r]

# Arrays for making system of equations
s = np.matlib.zeros(shape=(nx, 1))
c = np.matlib.zeros(shape=(nx, nx))

# Set up animation 
ylim = 1.1
fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(0, xmax*1e2), ylim=(0, ylim)) #units to cm
line, = ax.plot([], [], '-', color='#550000', lw=2)
time_template = '%.1f days'
time_text = ax.text(0.01, 0.96, '', ha='left', va='center', transform=plt.gca().transAxes, fontsize=14, weight='bold', color='black')
    
    
# System of equations solver of finite difference approximation of schrodinger equation
def animate(n):	
    if n == -1:
        line.set_data(x*1e2,r1*nx/sum_init) # units to cm
        time_text.set_text(time_template % (0))
    else:
        # Input Neumann Boundary Condition at x = 0
        c[0,0] = 1 + D*delta_t/(delta_x**2)
        s[0] = u[n,0]
        c[0,1] = -(D*delta_t)/(delta_x**2)
        
        # Input Open Boundary Condition at x = nx-1
        c[nx-1,nx-1] = 1
        s[nx-1] = u[n,nx-1]
        for j in xrange(1,nx-1):
            # Terms on RHS
            s[j] = u[n,j]
            # Coefficients on LHS with varying potential depending on position
            c[j,j-1] = -D*delta_t/(delta_x**2)
            c[j,j] = (1 + 2*D*delta_t/(delta_x**2))
            c[j,j+1] = -D*delta_t/(delta_x**2)
    
        # Compute the values of the solution using matrix division (i.e. cu = s -> u = c^-1*s)			
        sol = np.linalg.inv(c)*s
        # Input solutions into the solution array at time n+1 
        for k in xrange(0,nx):
            u[n+1, k] = sol[k]
        
        line.set_data(x*1e2, sol*nx/sum_init) # units to cm
        time_text.set_text(time_template % ((n+1)*delta_t))
    return line, time_text
        
ani = animation.FuncAnimation(fig, animate, np.arange(-1,nt-1),
                              interval=1, save_count=25, blit=True, init_func=init)
                           
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title(r'1D Diffusion of Molcules from a Meteorite',fontsize=14)
plt.ylabel('Fraction of total initial concentration',fontsize=14)
plt.xlabel(r'r (cm)',fontsize=14)

ani.save('1D_Diffusion.mp4', dpi = 300, writer=writer)
