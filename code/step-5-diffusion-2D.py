"""
 STEP FIVE: Diffusion in two dimensions.

 Our 2D diffusion equation is: 
 	du/dt = v * d^2u/dx^2 + v * d^2u/dy^2

 We discretize and solve for u^(n+1)_(i,j) to get:
 	u^(n+1)_(i,j) = u(^n)_(i,j) + ((v * del * t) / (del * x^2)) * (u(^n)_(i+1,j) - 2u(^n)_(i,j) + u^(n)_(i-1,j))
 					+ ((v * del * t) / (del * y^2)) * (u(^n)_(i,j+1) - 2u(^n)_(i,j) + u^(n)_(i-1,j-1))

 We still have the same initial conditions:
 	u, v = { 2 for x, y in (0.5, 1) x (0.5, 1)
 		   { 1 everywhere else

 And the same boundary conditions:
 	u = 1
 	v = 1 for x = 0, 2 and y = 0,2
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation, rc
from IPython.display import HTML
import time, sys
from mpl_toolkits.mplot3d import Axes3D

# Variable declarations
nx = 31
ny = 31
nt = 17
nu = .05
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .25
dt = sigma * dx * dy / nu

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)

u = np.ones((ny, nx)) 
un = np.ones((ny, nx))

# Assign initial conditions
u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2

# Plot initial conditions
fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')                      
X, Y = np.meshgrid(x, y)                            
surf = ax.plot_surface(X, Y, u, cmap=cm.plasma)

# Define a function for our diffusal phenomenon.
def diffuse(nt):
    u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2  
    
    for n in range(nt + 1): 
        un = u.copy()
        u[1:-1, 1:-1] = (un[1:-1, 1:-1] 
						+ nu * dt / dx**2 * (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2])
						+ nu * dt / dy**2 * (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]))
        u[0, :] = 1
        u[-1, :] = 1
        u[:, 0] = 1
        u[:, -1] = 1

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, u, rstride=1, cstride=1, cmap=cm.plasma,
        linewidth=0, antialiased=True)
    ax.set_zlim(1, 2.5)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$u$')
    ax.text2D(0.32, 0.95, '2D Diffusion at t=5', transform = ax.transAxes)
    plt.show()

diffuse(5)
