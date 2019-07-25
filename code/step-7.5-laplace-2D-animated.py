import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time, sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation, rc 
from IPython.display import HTML

# Variable declarations
nx = 31
ny = 31
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)

x = np.linspace(0, 2, nx)
y = np.linspace(0, 2, ny)

# Initial conditions
p = np.zeros((ny, nx))

# Boundary conditions
p[:, 0] = 0
p[:, -1] = y
p[0, :] = p[1, :]
p[-1, :] = p[-2, :]

# Initial figure
fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)
surf = ax.plot_surface(X, Y, p, rstride=1, cstride=1, cmap=cm.plasma, linewidth=0, antialiased=False)
ax.set_xlim(0, 2)
ax.set_ylim(0, 2)
ax.view_init(30, 225)
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_zlabel('$p$')
ax.text2D(0.30, .98, "2D Laplace Equation Over Time", transform = ax.transAxes)
plt.show()

# Initialization function for animation
def init():
	ax.clear()
	surf = ax.plot_surface(X, Y, p[:], rstride = 1, cstride = 1, cmap = cm.plasma, linewidth = 0, antialiased = False)
	return surf

# Animation function, similar to laplace2d function in step 7
def animate(i):
	pn = p.copy()
	p[1:-1, 1:-1] = ((dy**2 * (pn[1:-1, 2:] + pn[1:-1, 0:-2]) 
					+ dx**2 * (pn[2:, 1:-1] + pn[0:-2, 1:-1])) 
					/ (2 * (dx**2 + dy**2)))

	# Boundary conditions
	p[:, 0] = 0
	p[:, -1] = y
	p[0, :] = p[1, :]
	p[-1, :] = p[-2, :]

	ax.clear()
	surf = ax.plot_surface(X, Y, p[:], rstride = 1, cstride = 1, cmap = cm.plasma, linewidth = 0, antialiased = False)
	ax.set_xlim(0,2)
	ax.set_ylim(0,2)
	ax.view_init(30,225)
	ax.set_xlabel('$x$')
	ax.set_ylabel('$y$')
	ax.set_zlabel('$p$')
	ax.text2D(0.30, 0.98, "2D Laplace Equation Animated", transform = ax.transAxes)
	return surf

anim = animation.FuncAnimation(fig, animate, init_func = init, frames = 200, interval = 20)
anim.save('../LaplaceAnim.gif', writer = 'imagemagick', fps = 60)