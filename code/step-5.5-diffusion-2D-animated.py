"""
I thought it would be cool to see diffusion animated, so I adapted a lot of
the code from step 5 to make that happen. Instead of a diffuse function we
simply have an animation function, which does the same calculations and also
prints a figure. We can then call FuncAnimation to consecutively run our
animation function and give us our final gif. Overall, I'm happy with how
it turned out!
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation, rc 
from IPython.display import HTML
import time, sys
from mpl_toolkits.mplot3d import Axes3D

# Variable declarations
nx = 31
ny = 31
nt = 100
nu = .05
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .25
dt = sigma * dx * dy / nu

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)

u = np.ones((ny, nx)) 
un = np.ones((ny, nx))

# Initial conditions
u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

fig = plt.figure(figsize=(11, 7), dpi=100)
ax = fig.gca(projection='3d')                      
X, Y = np.meshgrid(x, y)                            
surf = ax.plot_surface(X, Y, u, cmap=cm.plasma)
ax.set_xlabel('$x$')
ax.set_zlabel('$u$')
ax.set_ylabel('$y$')
plt.show()

# Initialization function, called in FuncAnimation
def init():
	ax.clear()
	surf = ax.plot_surface(X, Y, u[:], cmap = cm.plasma)
	return surf

# Main animation function; similar to diffuse function in step 5, also prints figure
def animate(i):
    un = u.copy()
    u[1:-1, 1:-1] = (un[1:-1, 1:-1] 
					+ nu * dt / dx**2 * (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2])
					+ nu * dt / dy**2 * (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]))
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1

    ax.clear()
    surf = ax.plot_surface(X, Y, u[:], rstride=1, cstride=1, cmap=cm.plasma, linewidth= 0, antialiased=True)
    ax.set_zlim(1, 2.5)
    ax.set_xlabel('$x$')
    ax.set_zlabel('$u$')
    ax.set_ylabel('$y$')
    ax.text2D(0.35, 0.95, "2D Diffusion Over Time", transform=ax.transAxes);
    return surf

# Producing the animation
anim = animation.FuncAnimation(fig, animate, init_func = init, frames = nt, interval = 20)
anim.save('.../2Diff.gif', writer = 'imagemagick', fps = 60)
