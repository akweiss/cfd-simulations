"""
 STEP SIX: Burgers' equation in two dimensions.

 Our set of PDEs is as follows:
 	du/dt + u * du/dx + v * du/dy = v (d^2u/dx^2 +  d^2u/dy^2)
 	dv/dt + u * dv/dx + v * dv/dy = v (d^2v/dx^2 +  d^2v/dy^2)

 Here, we can copy our work from steps four and five to discretize the left 
 and right sides of our equation, respectively. Doing this and solving for
 u^(n+1)_(i,j) and v^(n+1)_(i,j) we get:
  	u^(n+1)_(i,j) = u(^n)_(i,j) - u(^n)_(i,j) * ((del * t) / (del * x)) * (u(^n)_(i,j) - u(^n)_(i-1,j)) 
  					- v(^n)_(i,j) * ((del * t) / (del * y)) * (u(^n)_(i,j) - u(^n)_(i,j-1)) 
  					+ ((v * del * t) / (del * x^2)) * (u(^n)_(i+1,j) - 2u(^n)_(i,j) + u^(n)_(i-1,j)) 
  					+ ((v * del * t) / (del * y^2)) * (u(^n)_(i,j+1) - 2u(^n)_(i,j) + u^(n)_(i-1,j-1))

  	v^(n+1)_(i,j) = v(^n)_(i,j) - u(^n)_(i,j) * ((del * t) / (del * x)) * (v(^n)_(i,j) - v(^n)_(i-1,j)) 
  					- v(^n)_(i,j) * ((del * t) / (del * y)) * (v(^n)_(i,j) - v(^n)_(i,j-1)) 
  					+ ((v * del * t) / (del * x^2)) * (v(^n)_(i+1,j) - 2v(^n)_(i,j) + v^(n)_(i-1,j)) 
  					+ ((v * del * t) / (del * y^2)) * (v(^n)_(i,j+1) - 2v(^n)_(i,j) + v^(n)_(i-1,j-1))

 Here, our initial conditions are:
 	u, v = { 2 for x, y in (0.5, 1) x (0.5, 1)
 			{ 1 everywhere else

 And our boundary conditions are:
 	u = 1
 	v = 1 for x = 0, 2 and y = 0,2
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time, sys
from mpl_toolkits.mplot3d import Axes3D

# Variable declarations
nx = 41
ny = 41
nt = 120
nu = .01
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .0009
dt = sigma * dx * dy / nu

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)

u = np.ones((ny, nx)) 
un = np.ones((ny, nx))
v = np.ones((ny, nx)) 
vn = np.ones((ny, nx)) 

# Assign initial conditions
u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2
v[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2

# Plot initial conditions
fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')                      
X, Y = np.meshgrid(x, y)                            
surf = ax.plot_surface(X, Y, u, cmap=cm.plasma)

# Implement our solutions for u^(n+1)_(i,j) and v^(n+1)_(i,j) using array operations.
for n in range(nt + 1): 
    un = u.copy()
    vn = v.copy()
 
    u[1:-1, 1:-1] = (un[1:-1, 1:-1] 
					- dt / dx * un[1:-1, 1:-1] * (un[1:-1, 1:-1] - un[1:-1, 0:-2])
					- dt / dy * vn[1:-1, 1:-1] * (un[1:-1, 1:-1] - un[0:-2, 1:-1])
					+ nu * dt / dx**2 * (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2])
					+ nu * dt / dy**2 * (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]))

    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] 
					- dt / dx * un[1:-1, 1:-1] * (vn[1:-1, 1:-1] - vn[1:-1, 0:-2])
					- dt / dy * vn[1:-1, 1:-1] * (vn[1:-1, 1:-1] - vn[0:-2, 1:-1])
					+ nu * dt / dx**2 * (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2])
					+ nu * dt / dy**2 * (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1]))

    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1

    v[0, :] = 1
    v[-1, :] = 1
    v[:, 0] = 1
    v[:, -1] = 1

fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, u, rstride=1, cstride=1, cmap=cm.plasma, linewidth=0, antialiased=True)
ax.set_zlim(1, 2.5)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$');
plt.show()
