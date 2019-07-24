"""
 STEP NINE: Cavity flow with Navier-Stokes.

 We have 3 PDEs in this system: two equations for velocity, and one equation
 for pressure. They are as follows:
 	du/dt + u * du/dx + v * du/dy = - (1 / rho) dp/dx + v * (d^2u/dx^2 + d^2u/dy^2)
 	dv/dt + u * dv/dx + v * dv/dy = - (1 / rho) dp/dx + v * (d^2v/dx^2 + d^2v/dy^2)
 	d^2p/dx^2 + d^2p/dy^2 = - rho * ((du/dx * du/dx) + (2 * du/dy * dv/dx) + (dv/dy * dv/dy))

 We then discretize and solve. Doing this, we get the momentum equation in the u direction:
 	u(^n+1)_(i,j) = u(^n)_(i,j) - u(^n)_(i,j) * ((del * t) / (del * x)) * ((u(^n)_(i,j) - u(^n)_(i-1,j)))
 					- v(^n)_(i,j) * ((del * t) / (del * y)) * ((u(^n)_(i,j) - u(^n)_(i,j-1)))
 					- ((del * t) / (rho * 2 * del * x)) * ((p(^n)_(i+1,j) - p(^n)_(i-1,j)))
 					+ v * [((del * t) / (del * x^2)) * ((u(^n)_(i+1,j) - 2u(^n)_(i,j) + u(^n)_(i-1,j)))
 					((del * t) / (del * y^2)) * ((u(^n)_(i,j+1) - 2u(^n)_(i,j) + u(^n)_(i,j-1)))]

 The momentum equation in the v direction:
 	v(^n+1)_(i,j) = v(^n)_(i,j) - u(^n)_(i,j) * ((del * t) / (del * x)) * ((v(^n)_(i,j) - v(^n)_(i-1,j)))
 					- v(^n)_(i,j) * ((del * t) / (del * y)) * ((v(^n)_(i,j) - v(^n)_(i,j-1)))
 					- ((del * t) / (rho * 2 * del * x)) * ((p(^n)_(i+1,j) - p(^n)_(i-1,j)))
 					+ v * [((del * t) / (del * x^2)) * ((v(^n)_(i+1,j) - 2v(^n)_(i,j) + v(^n)_(i-1,j)))
 					((del * t) / (del * y^2)) * ((v(^n)_(i,j+1) - 2v(^n)_(i,j) + v(^n)_(i,j-1)))]

 And the pressure (Poisson) equation:
 	p(^n)_(i,j) = [(p(^n)_(i+1,j) + p(^n)_(i-1,j) * del * y^2) + (p(^n)_(i,j+1) + p(^n)_(i,j+1) * del * x^2)] 
 				  / (2 (del * x^2 + del * y^2))
 				  - ((rho * del * x^2 * del * y^2)) / (2 (del * x^2 + del * y^2))
 				  * [ ((1 / del * t) * ((u_(i+1,j) - u_(i-1,j) / (2 * del * x)) + ((v_(i,j+1) - v_(i,j-1) / (2 * del * y)))
 				  - ((u_(i+1,j) - u_(i-1,j) / (2 * del * x)) * ((u_(i+1,j) - u_(i-1,j) / (2 * del * x))
 				  - 2 * ((u_(i,j+1) - u_(i,j-1) / (2 * del * y)) * ((v_(i+1,j) - v_(i-1,j) / (2 * del * x))
 				  - ((v_(i,j+1) - v_(i,j-1) / (2 * del * y)) * ((v_(i,j-1) - v_(i,j-1) / (2 * del * y)) ]

 Here, our initial conditions are:
 	u, v, p = 0 everywhere

 And our boundary conditions are:
 	u = 1 at y = 2
 	u, v = 0 on other boundaries
 	dp/dy = 0 at y = 0
 	p = 0 at y = 2
 	dp/dx = 0 at x = 0, 2

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time, sys
from mpl_toolkits.mplot3d import Axes3D

# Variable declarations
nx = 41
ny = 41
nt = 500
nit = 50
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)

x = np.linspace(0, 2, nx)
y = np.linspace(0, 1, ny)
X, Y = np.meshgrid(x, y)

rho = 1
nu = .1
dt = .001

u = np.zeros((ny,nx))
v = np.zeros((ny,nx))
p = np.zeros((ny,nx))
b = np.zeros((ny,nx))

# We'll use two functions to break up the calculations of the pressure Poisson equation, for the sake of clarity and managing errors.

# The first function deals with the second half of the equation, which I bracketed in my overview
def build_up(b, rho, dt, u, v, dx, dy):
	b[1:-1, 1:-1] = (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) 
					- ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 
					- 2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) * (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) 
					- ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))
	return b

# The second function handles the remainder of the Poisson equation and boundary conditions
def poisson(p, dx, dy, b):
	pn = np.empty_like(p)
	pn = p.copy()

	for q in range(nit):
		pn = p.copy()
		p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 + (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) 
						/ (2 * (dx**2 + dy**2)) 
						- dx**2 * dy**2 / (2 * (dx**2 + dy**2)) 
						* b[1:-1,1:-1])
		# Boundary conditions
		p[:, -1] = p[:, -2]			# dp/dx = 0 at x = 2
		p[0, :] = p[1, :]			# dp/dy = 0 at y = 0
		p[:, 0] = p[:, 1]			# dp/dx = 0 at x = 0
		p[-1, :] = 0				# p = 0 at y = 2

	return p



# This function handles both momentum equations and their boundary conditions, 
# incorporating our p from previous Poisson functions for convenient plotting
def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu):
	un = np.empty_like(u)
	vn = np.empty_like(v)
	b = np.zeros((ny,nx))

	for n in range(nt):
		un = u.copy()
		vn = v.copy()

		b = build_up(b, rho, dt, u, v, dx, dy)
		p = poisson(p, dx, dy, b)

		u[1:-1, 1:-1] = (un[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx * (un[1:-1, 1:-1] - un[1:-1, 0:-2]) 
						- vn[1:-1, 1:-1] * dt / dy * (un[1:-1, 1:-1] - un[0:-2, 1:-1])
						- dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) 
						+ nu * (dt / dx**2 * (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) 
						+ dt / dy**2 * (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))
		v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx * (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) 
						- vn[1:-1, 1:-1] * dt / dy * (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) 
						- dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) 
						+ nu * (dt / dx**2 * (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) 
						+ dt / dy**2 * (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

		# Boundary conditions
		u[0, :] = 0
		u[:, 0] = 0
		u[:, -1] = 0
		u[-1, :] = 1		# Velocity on cavity lid is 1
		v[0, :] = 0
		v[:, 0] = 0
		v[:, -1] = 0
		v[-1, :] = 0

	return u, v, p

# Conditions for the solver
u = np.zeros((ny,nx))
v = np.zeros((ny,nx))
p = np.zeros((ny,nx))
b = np.zeros((ny,nx))
nt = 700
u, v, p = cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu)

# Plot the results as a velocity field
fig = plt.figure(figsize=(11,7), dpi=100)
plt.contourf(X, Y, p, alpha=0.5, cmap=cm.plasma)  
plt.colorbar()
plt.contour(X, Y, p, cmap=cm.plasma)
plt.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2]) 
plt.xlabel('X')
plt.ylabel('Y');

# We can also visualize with a streamplot
fig = plt.figure(figsize=(11, 7), dpi=100)
plt.contourf(X, Y, p, alpha=0.5, cmap=cm.plasma)
plt.colorbar()
plt.contour(X, Y, p, cmap=cm.plasma)
plt.streamplot(X, Y, u, v)
plt.xlabel('X')
plt.ylabel('Y');
plt.show()
