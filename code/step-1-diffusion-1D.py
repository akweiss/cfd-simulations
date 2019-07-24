"""
STEP ONE: Diffusion equation in one dimension.

The diffusion equation in one dimension is:
	du/dt = v * d^2u/dx^2

We need to discretize this second-order derivative using a central difference scheme, 
then solve for the unknown, u(^n+1)_i. In this case,
	u(^n+1)_(i) = u(^n)_i + ((v * del * t) / (del * x^2)) * (u(^n)_(i+1) - 2u(^n)_i + u(^n)_(i-1))

Then we apply this to our solver.
"""

import numpy
import matplotlib.pyplot as plt

nx = 41						# Discretization of the domain
dx = 2 / (nx - 1)
nt = 20    					# Number of timesteps 
nu = 0.3   					# Viscosity
sigma = .2 					
dt = sigma * dx**2 / nu 	# Amount of time each timestep covers

u = numpy.ones(nx)      # All elements are 1
u[int(.5 / dx):int(1 / dx + 1)] = 2  # u = 2 between 0.5 and 1, per initial conditions

un = numpy.ones(nx) # Temporary array

for n in range(nt):  	# Iterates nt times
    un = u.copy() 		# Copies values of u into un
    for i in range(1, nx - 1):
    	# 2nd derivative finite difference application
        u[i] = un[i] + nu * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1])		

plt.plot(numpy.linspace(0, 2, nx), u)
plt.show()
