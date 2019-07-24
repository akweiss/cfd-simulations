"""
 STEP FOUR: Nonlinear convection in two dimensions.

 In this step we want to solve 2D nonlinear convection, which is represented
 by these corresponding PDEs:
 	du/dt + u * du/dx + v * du/dy = 0
 	dv/dt + u * dv/dx + v * dv/dy = 0

 We can discretize and solve these in a similar fashion to step three and get:
 	u^(n+1)_(i,j) = u(^n)_(i,j) - u_(i,j) * ((del * t) / (del * x)) * (u(^n)_(i,j) - u(^n)_(i-1,j)) 
 					- v(^n)_(i,j) * ((del * t) / (del * y)) * (u(^n)_(i,j) - u(^n)_(i,j-1))

 	v^(n+1)_(i,j) = v(^n)_(i,j) - u_(i,j) * ((del * t) / (del * x)) * (u(^n)_(i,j) - u(^n)_(i-1,j)) 
 					- v(^n)_(i,j) * ((del * t) / (del * y)) * (u(^n)_(i,j) - u(^n)_(i,j-1))

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
nx = 101
ny = 101
nt = 100
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .2
dt = sigma * dx

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
    u[1:, 1:] = (un[1:, 1:] - (un[1:, 1:] * c * dt / dx * (un[1:, 1:] - un[1:, :-1]))
    						- vn[1:, 1:] * c * dt / dy * (un[1:, 1:] - un[:-1, 1:]))

    v[1:, 1:] = (vn[1:, 1:] - (un[1:, 1:] * c * dt / dx * (vn[1:, 1:] - vn[1:, :-1]))
    						- vn[1:, 1:] * c * dt / dy * (vn[1:, 1:] - vn[:-1, 1:]))

    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1

    v[0, :] = 1
    v[-1, :] = 1
    v[:, 0] = 1
    v[:, -1] = 1

# Plot the PDEs
fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf2 = ax.plot_surface(X, Y, u, cmap=cm.plasma)
ax.set_xlabel('$x$')
ax.set_zlabel('$u$')
ax.set_ylabel('$y$')
ax.text2D(0.35, 0.95, "2D Non-Linear Convection at t=10 for u Velocity", transform=ax.transAxes)

fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf3 = ax.plot_surface(X, Y, v, cmap=cm.plasma)
ax.set_xlabel('$x$')
ax.set_zlabel('$v$')
ax.set_ylabel('$y$')
ax.text2D(0.35, 0.95, "2D Non-Linear Convection at t=10 for v Velocity", transform=ax.transAxes)

plt.show()
