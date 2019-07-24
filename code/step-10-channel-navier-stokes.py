"""
 STEP TEN: Channel flow with Navier-Stokes.

 This is nearly identical to step nine, but we are adding a source term, F, for the u momentum
 equation. Doing so, our system of Navier-Stokes equations therefore becomes:
  	du/dt + u * du/dx + v * du/dy = - (1 / rho) dp/dx + v * (d^2u/dx^2 + d^2u/dy^2) + F
 	dv/dt + u * dv/dx + v * dv/dy = - (1 / rho) dp/dx + v * (d^2v/dx^2 + d^2v/dy^2)
 	d^2p/dx^2 + d^2p/dy^2 = - rho * ((du/dx * du/dx) + (2 * du/dy * dv/dx) + (dv/dy * dv/dy))

 Now we discretize and solve both velocity equations at time n+1.

 The momentum equation in the u direction:
  	u(^n+1)_(i,j) = u(^n)_(i,j) - u(^n)_(i,j) * ((del * t) / (del * x)) * ((u(^n)_(i,j) - u(^n)_(i-1,j)))
 					- v(^n)_(i,j) * ((del * t) / (del * y)) * ((u(^n)_(i,j) - u(^n)_(i,j-1)))
 					- ((del * t) / (rho * 2 * del * x)) * ((p(^n)_(i+1,j) - p(^n)_(i-1,j)))
 					+ v * [((del * t) / (del * x^2)) * ((u(^n)_(i+1,j) - 2u(^n)_(i,j) + u(^n)_(i-1,j)))
 					((del * t) / (del * y^2)) * ((u(^n)_(i,j+1) - 2u(^n)_(i,j) + u(^n)_(i,j-1)))]
 					+ del * t * F

 The momentum equation in the v direction:
  	v(^n+1)_(i,j) = v(^n)_(i,j) - u(^n)_(i,j) * ((del * t) / (del * x)) * ((v(^n)_(i,j) - v(^n)_(i-1,j)))
 					- v(^n)_(i,j) * ((del * t) / (del * y)) * ((v(^n)_(i,j) - v(^n)_(i,j-1)))
 					- ((del * t) / (rho * 2 * del * x)) * ((p(^n)_(i+1,j) - p(^n)_(i-1,j)))
 					+ v * [((del * t) / (del * x^2)) * ((v(^n)_(i+1,j) - 2v(^n)_(i,j) + v(^n)_(i-1,j)))
 					((del * t) / (del * y^2)) * ((v(^n)_(i,j+1) - 2v(^n)_(i,j) + v(^n)_(i,j-1)))]

 And the pressure equation:
  	p(^n)_(i,j) = [(p(^n)_(i+1,j) + p(^n)_(i-1,j) * del * y^2) + (p(^n)_(i,j+1) + p(^n)_(i,j+1) * del * x^2)] 
 				  / (2 (del * x^2 + del * y^2))
 				  - ((rho * del * x^2 * del * y^2)) / (2 (del * x^2 + del * y^2))
 				  * [ ((1 / del *t) * ((u_(i+1,j) - u_(i-1,j) / (2 * del * x)) + ((v_(i,j+1) - v_(i,j-1) / (2 * del * y)))
 				  - ((u_(i+1,j) - u_(i-1,j) / (2 * del * x)) * ((u_(i+1,j) - u_(i-1,j) / (2 * del * x))
 				  - 2 * ((u_(i,j+1) - u_(i,j-1) / (2 * del * y)) * ((v_(i+1,j) - v_(i-1,j) / (2 * del * x))
 				  - ((v_(i,j+1) - v_(i,j-1) / (2 * del * y)) * ((v_(i,j-1) - v_(i,j-1) / (2 * del * y)) ]

 Note these are identical to the equations from step nine, except the momentum equation in the u direction, 
 which now has the additional del * t * F term.

 Our initial conditions are:
 	u, v, p = 0 everywhere

 And our boundary conditions are:
 	u, v, p are periodic on x = 0, 2
 	u, v = 0 at y = 0, 2
 	dp/dy = 0 at y = 0, 2
 	F = 1 everywhere
 
 Therefore, our main task is to adapt the code from step nine to accomodate these periodic boundary conditions.
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
F = 1

u = np.zeros((ny,nx))
v = np.zeros((ny,nx))
p = np.zeros((ny,nx))
b = np.zeros((ny,nx))

# We'll use two functions to break up the calculations of the pressure Poisson equation, for the sake of clarity and managing errors.

# The first function deals with the second half of the equation, accounting for the periodic boundary conditions at x = 0, 2

def build_up(rho, dt, dx, dy, u, v):
	b = np.zeros_like(u)
	b[1:-1, 1:-1] = (rho * (1 / dt * ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) 
					- ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx))**2 
					- 2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) * (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx))
					- ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy))**2))
    
    # Periodic boundary condition at x = 2
	b[1:-1, -1] = (rho * (1 / dt * ((u[1:-1, 0] - u[1:-1,-2]) / (2 * dx) + (v[2:, -1] - v[0:-2, -1]) / (2 * dy)) 
				- ((u[1:-1, 0] - u[1:-1, -2]) / (2 * dx))**2 
				- 2 * ((u[2:, -1] - u[0:-2, -1]) / (2 * dy) * (v[1:-1, 0] - v[1:-1, -2]) / (2 * dx)) 
				- ((v[2:, -1] - v[0:-2, -1]) / (2 * dy))**2))

    # Periodic boundary condition at x = 0
	b[1:-1, 0] = (rho * (1 / dt * ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx) + (v[2:, 0] - v[0:-2, 0]) / (2 * dy)) 
				- ((u[1:-1, 1] - u[1:-1, -1]) / (2 * dx))**2 
				- 2 * ((u[2:, 0] - u[0:-2, 0]) / (2 * dy) * (v[1:-1, 1] - v[1:-1, -1]) / (2 * dx))
				- ((v[2:, 0] - v[0:-2, 0]) / (2 * dy))**2))
	return b

# The second function handles the remainder of the Poisson equation, the periodic boundary conditions, and wall boundaries

def poisson(p, dx, dy):
    pn = np.empty_like(p)
    
    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy**2 + (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx**2) 
        				/ (2 * (dx**2 + dy**2)) 
        				- dx**2 * dy**2 / (2 * (dx**2 + dy**2)) 
        				* b[1:-1, 1:-1])

        # Periodic boundary condition at x = 2
        p[1:-1, -1] = (((pn[1:-1, 0] + pn[1:-1, -2])* dy**2 + (pn[2:, -1] + pn[0:-2, -1]) * dx**2) 
        			/ (2 * (dx**2 + dy**2)) 
        			- dx**2 * dy**2 / (2 * (dx**2 + dy**2)) 
        			* b[1:-1, -1])

        # Periodic boundary condition x = 0
        p[1:-1, 0] = (((pn[1:-1, 1] + pn[1:-1, -1])* dy**2 + (pn[2:, 0] + pn[0:-2, 0]) * dx**2) 
        			/ (2 * (dx**2 + dy**2)) 
        			- dx**2 * dy**2 / (2 * (dx**2 + dy**2)) 
        			* b[1:-1, 0])
        
        # Wall boundary conditions
        p[-1, :] =p[-2, :]	# dp/dy = 0 at y = 2
        p[0, :] = p[1, :]	# dp/dy = 0 at y = 0
    
    return p


# The following handles both momentum equations and their boundary conditions,
# using an approach similar to step 7 where we iterate until the difference 
# between two consecutive steps is very small.

udiff = 1
stepcount = 0

while udiff > .001:
	un = u.copy()
	vn = v.copy()

	b = build_up(rho, dt, dx, dy, u, v)
	p = poisson(p, dx, dy)

	u[1:-1, 1:-1] = (un[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx * (un[1:-1, 1:-1] - un[1:-1, 0:-2]) 
					- vn[1:-1, 1:-1] * dt / dy * (un[1:-1, 1:-1] - un[0:-2, 1:-1]) 
					- dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) 
					+ nu * (dt / dx**2 * (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) 
					+ dt / dy**2 * (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])) + F * dt)

	v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - un[1:-1, 1:-1] * dt / dx * (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) 
					- vn[1:-1, 1:-1] * dt / dy * (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) 
					- dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) 
					+ nu * (dt / dx**2 * (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) 
					+ dt / dy**2 * (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

    # Periodic boundary condition u at x = 2     
	u[1:-1, -1] = (un[1:-1, -1] - un[1:-1, -1] * dt / dx * (un[1:-1, -1] - un[1:-1, -2]) 
				- vn[1:-1, -1] * dt / dy * (un[1:-1, -1] - un[0:-2, -1]) 
				- dt / (2 * rho * dx) * (p[1:-1, 0] - p[1:-1, -2]) 
				+ nu * (dt / dx**2 * (un[1:-1, 0] - 2 * un[1:-1,-1] + un[1:-1, -2]) 
				+ dt / dy**2 * (un[2:, -1] - 2 * un[1:-1, -1] + un[0:-2, -1])) + F * dt)

    # Periodic boundary condition u at x = 0
	u[1:-1, 0] = (un[1:-1, 0] - un[1:-1, 0] * dt / dx * (un[1:-1, 0] - un[1:-1, -1]) 
				- vn[1:-1, 0] * dt / dy * (un[1:-1, 0] - un[0:-2, 0]) 
				- dt / (2 * rho * dx) * (p[1:-1, 1] - p[1:-1, -1]) 
				+ nu * (dt / dx**2 * (un[1:-1, 1] - 2 * un[1:-1, 0] + un[1:-1, -1]) 
				+ dt / dy**2 * (un[2:, 0] - 2 * un[1:-1, 0] + un[0:-2, 0])) + F * dt)

    # Periodic boundary condition v at x = 2
	v[1:-1, -1] = (vn[1:-1, -1] - un[1:-1, -1] * dt / dx * (vn[1:-1, -1] - vn[1:-1, -2]) 
				- vn[1:-1, -1] * dt / dy * (vn[1:-1, -1] - vn[0:-2, -1]) 
				- dt / (2 * rho * dy) * (p[2:, -1] - p[0:-2, -1]) 
				+ nu * (dt / dx**2 * (vn[1:-1, 0] - 2 * vn[1:-1, -1] + vn[1:-1, -2]) 
				+ dt / dy**2 * (vn[2:, -1] - 2 * vn[1:-1, -1] + vn[0:-2, -1])))

    # Periodic boundary condition v at x = 0
	v[1:-1, 0] = (vn[1:-1, 0] - un[1:-1, 0] * dt / dx * (vn[1:-1, 0] - vn[1:-1, -1]) 
				- vn[1:-1, 0] * dt / dy * (vn[1:-1, 0] - vn[0:-2, 0]) 
				- dt / (2 * rho * dy) * (p[2:, 0] - p[0:-2, 0]) 
				+ nu * (dt / dx**2 * (vn[1:-1, 1] - 2 * vn[1:-1, 0] + vn[1:-1, -1]) 
				+ dt / dy**2 * (vn[2:, 0] - 2 * vn[1:-1, 0] + vn[0:-2, 0])))


    # Wall boundary conditions: u, v = 0 at y = 0, 2
	u[0, :] = 0
	u[-1, :] = 0
	v[0, :] = 0
	v[-1, :]=0
    
	udiff = (np.sum(u) - np.sum(un)) / np.sum(u)
	stepcount += 1


# Plot the results
fig = plt.figure(figsize = (11,7), dpi=100)
plt.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Channel Flow with Periodic Boundary Conditions')
plt.show()
