"""
 STEP EIGHT: Poisson equation in two dimensions.

 When we have incompressible flow, a kinematic constraint arises in
 the equation that represents mass conservation at constant density.
 The constraint requires the pressure field to evolve such that the
 rate of expansion vanishes everywhere. To avoid this, we can instead
 construct our own pressure field that guarantees continuity is sati-
 sfied by taking the divergence of the momentum equation. Doing this,
 we arrive at the Poisson equation for pressure.

 Our 2D Poisson equation with b as the source term is: 
 	b = d^2p/dx^2 + d^2p/dy^2

 We can discretize this like we have in the previous steps to get:
 	p(^n)_(i,j) = [del * y^2 (p(^n)_(i+1,j) + p(^n)_(i-1,j)) 
 				   + del * x^2 (p(^n)_(i,j+1) + p(^n)_(i,j-1))
 				   - (del * x^2) * (del * y^2) * (b(^n)_(i,j))] 
 				   / [2 (del * x^2 + del * y^2)]

 We will solve this assuming p = 0 everywhere initially, and applying
 the following boundary conditions:
 	p = 0 at x = 0, 2 and y = 0, 1

 The source term has two initial spikes in the domain, as follows:
 	b_(i,j) = 100 at i = (1/4) * nx, j = (1/4) * ny
 	b_(i,j) = -100 at i = (3/4) * nx, j = (3/4) * ny
 	b_(i,j) = 0 everywhere else

"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time, sys
from mpl_toolkits.mplot3d import Axes3D

# Variable declarations
nx = 50
ny = 50
nt  = 100
xmin = 0
xmax = 2
ymin = 0
ymax = 1

dx = (xmax - xmin) / (nx - 1)
dy = (ymax - ymin) / (ny - 1)

# Initialization
p  = np.zeros((ny, nx))
pd = np.zeros((ny, nx))
b  = np.zeros((ny, nx))
x  = np.linspace(xmin, xmax, nx)
y  = np.linspace(xmin, xmax, ny)

# Source term initial spikes
b[int(ny / 4), int(nx / 4)]  = 100
b[int(3 * ny / 4), int(3 * nx / 4)] = -100

#crazy shit here
for it in range(nt):
	pd = p.copy()
	p[1:-1, 1:-1] = (((pd[1:-1, 2:] + pd[1:-1, :-2]) * dy**2 
					+ (pd[2:, 1:-1] + pd[:-2, 1:-1]) * dx**2 
					- b[1:-1, 1:-1] * dx**2 * dy**2) 
					/ (2 * (dx**2 + dy**2)))
	p[0, :] = 0
	p[ny-1, :] = 0
	p[:, 0] = 0
	p[:, nx-1] = 0

# We can reuse our plotting function from step 7 to plot our results
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

plot2d(x, y, p)

