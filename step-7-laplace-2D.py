"""
 STEP SEVEN: Laplace's equation in two dimensions.

 Our equation is as follows: 
 	d^2p/dx^2 + d^2p/dy^2 = 0

 Because Laplace's equation has features reminiscent of diffusion, we should discretize our
 PDE using central differences - as opposed to backward or forward diffences - in order to
 accurately simulate the physics of this phenomenon. Additionally, the Laplace equation does
 not have time dependence; it calculates the equilibrium state of a system given a set of
 boundary conditions.

 Because of this, our approach from steps 5-8 alters somewhat. Time independence means that
 instead of calculating where some system will be at time t, we must solve for p(^n)_(i,j)
 until it meets some condition that we specify. In this case, we can approximate an equili-
 brium state when the difference between two consecutive iterative steps is very small.

 Discretizing using central difference and solving for p(^n)_(i,j) we get:
 	p(^n)_(i,j) = [del * y^2 (p(^n)_(i+1,j) + p(^n)+_(i-1,j)) 
 				  + del * x^2 (p(^n)_(i,j+1) + p(^n)+_(i,j-1))] 
 				  / [2 (del * x^2 + del * y^2)]

 We will solve the Laplace equation numerically, assuming initially that p = 0 everywhere. 
 Then we will employ the following boundary conditions:
 	p = 0 at x = 0,
 	p = y at x = 2,
 	dp/dy = 0 at y = 0, 1
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time, sys
from mpl_toolkits.mplot3d import Axes3D

# We'll define two functions for this step. The first we'll use to plot our results.
def plot2d(x, y, p):
	fig = plt.figure()
	ax = fig.gca(projection='3d')
	X, Y = np.meshgrid(x, y)
	surf = ax.plot_surface(X, Y, p, rstride=1, cstride=1, cmap=cm.plasma, linewidth=0, antialiased=True)
	ax.set_xlim(0, 2)
	ax.set_ylim(0, 2)
	ax.view_init(30, 225)
	ax.set_xlabel('$x$')
	ax.set_ylabel('$y$');
	plt.show()


# Our second function we'll use to solve our PDE.
def laplace2d(p, y, dx, dy, target):
	norm = 1
	pn = np.empty_like(p)

	while norm > target:
		pn = p.copy()
		# Solution for p(^n)_(i,j)
		p[1:-1, 1:-1] = ((dy**2 * (pn[1:-1, 2:] + pn[1:-1, 0:-2]) 
						+ dx**2 * (pn[2:, 1:-1] + pn[0:-2, 1:-1])) 
						/ (2 * (dx**2 + dy**2)))
		# Boundary conditions
		p[:, 0] = 0
		p[:, -1] = y
		p[0, :] = p[1, :]
		p[-1, :] = p[-2, :]
		norm = (np.sum(np.abs(p[:]) - np.abs(pn[:])) / np.sum(np.abs(pn[:])))

	return p

# Variable declarations
nx = 31
ny = 31
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)

x = np.linspace(0, 2, nx)
y = np.linspace(0, 1, ny)

# Initial conditions
p = np.zeros((ny, nx))

# Boundary conditions
p[:, 0] = 0
p[:, -1] = y
p[0, :] = p[1, :]
p[-1, :] = p[-2, :]

# Plot initial conditions using our plotting function
plot2d(x, y, p)

# Run Laplace function, with a target of .01
p = laplace2d(p, y, dx, dy, 1e-4)

# Plot the new value of p using our plotting function
plot2d(x, y, p)