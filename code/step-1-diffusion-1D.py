# STEP ONE: Diffusion equation in one dimension.

import numpy
import matplotlib.pyplot as plt

nx = 41						# Discretization of the domain
dx = 2 / (nx - 1)
nt = 20    					# Number of timesteps 
nu = 0.3   					# Viscosity
sigma = .2 					
dt = sigma * dx**2 / nu 	# Amount of time each timestep covers

# Auxilary variables for difference formulas
for i in 1:nx:
	ip1(i) = i+1
	im1(i) = i-1
	x(i) = (i-1)*dx
ip1(nx) = 1
im1(1) = nx


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
