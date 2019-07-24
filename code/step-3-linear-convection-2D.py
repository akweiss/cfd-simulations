"""
 STEP THREE: Linear convection in two dimensions.

 We can expand the examples of linear convection, nonlinear convection, diffusion,
 and Burgers' equation into two dimensions by directly applying the definition of
 partial derivatives. That is, a partial derivative with respect to x is the vari-
 ation in the x direction at constant y.

Linear convection in 2D can be written as:
 	du/dt + c * du/dx + c * du/dy = 0

 We will discretize the timestep dt with forward difference approximations, and the
 two spatial steps will be discretized using backward difference approximations.

 Thus, we can discretize our PDE and solve for our only unknown, u(^(n+1))_(i,j).
 	u^(n+1)_(i,j) = u(^n)_(i,j) - c * ((del * t) / (del * x)) * (u(^n)_(i,j) - u(^n)_(i-1,j))
 					- c * ((del * t) / (del * y)) * (u(^n)_(i,j) - u(^n)_(i,j-1))
 This is the main component we will be implementing.

 For this equation, we will also use the following initial conditions:
 	u(x,y) = { 2 for 0.5 <= x,y <= 2
 			 { 1 for everywhere else 

 And the following boundary conditions:
 	u = 1 for x = 0, 2 and y = 0, 2
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

# Assign initial conditions
# Recall that our initial conditions are:
# 	u(x,y) = 2 for 0.5 <= x,y <= 2
#			 1 everywhere else
u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2 

# Plot initial conditions
fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u, cmap=cm.plasma)

# We can use array operations here to eliminate having to use nested for loops
# in order to evaluate our wave in two dimensions. 
for n in range(nt + 1): 
    un = u.copy()
    # Here is where we implement our result for u(^n+1)_(i,j) after discretizing and solving.
    u[1:, 1:] = (un[1:, 1:] - (c * dt / dx * (un[1:, 1:] - un[1:, :-1])) -
                              (c * dt / dy * (un[1:, 1:] - un[:-1, 1:])))
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1

# Plot the PDE
fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')
surf2 = ax.plot_surface(X, Y, u, cmap=cm.plasma)
ax.set_zlabel('$u$')
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.text2D(0.3, 0.95, 'Linear Convection in 2D at t=10', transform=ax.transAxes)
plt.show()
