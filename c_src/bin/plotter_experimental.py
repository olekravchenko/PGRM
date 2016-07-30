#!/usr/bin/python
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np
import sys

fig = plt.figure(figsize=plt.figaspect(0.7))
#ax = fig.gca(projection='3d')

##First Sub-Plot ##

ax = fig.add_subplot(2, 2, 1, projection='3d')
plt.title('Numeric Solution')
data = np.genfromtxt('../plot_data/plot_region')
x = data[:,0]
y = data[:,1]
z = data[:,2]

xi = np.linspace(min(x), max(x))
yi = np.linspace(min(y), max(y))

X, Y = np.meshgrid(xi, yi)
Z = griddata(x, y, z, xi, yi)

surf = ax.plot_surface(X, Y, Z, rstride=2, cstride=2, cmap=cm.jet,
                       linewidth=0.02, antialiased=True)

ax.set_zlim3d(np.min(Z), np.max(Z))
fig.colorbar(surf)

##Second Sub-Plot ##

ax = fig.add_subplot(2, 2, 2, projection='3d')
plt.title('Error')
data = np.genfromtxt('../plot_data/plot_region_error')
x = data[:,0]
y = data[:,1]
z = data[:,2]

xi = np.linspace(min(x), max(x))
yi = np.linspace(min(y), max(y))

X, Y = np.meshgrid(xi, yi)
Z = griddata(x, y, z, xi, yi)

surf = ax.plot_surface(X, Y, Z, rstride=2, cstride=2, cmap=cm.jet,
                       linewidth=0.02, antialiased=True)

ax.set_zlim3d(np.min(Z), np.max(Z))
fig.colorbar(surf)

##Third Sub-Plot ##

ax = fig.add_subplot(2, 2, 3, projection='3d')
plt.title('Exact Solution')
data = np.genfromtxt('../plot_data/plot_exact_solution')
x = data[:,0]
y = data[:,1]
z = data[:,2]

xi = np.linspace(min(x), max(x))
yi = np.linspace(min(y), max(y))

X, Y = np.meshgrid(xi, yi)
Z = griddata(x, y, z, xi, yi)

surf = ax.plot_surface(X, Y, Z, rstride=2, cstride=2, cmap=cm.jet,
                       linewidth=0.02, antialiased=True)

ax.set_zlim3d(np.min(Z), np.max(Z))
fig.colorbar(surf)

##Second Sub-Plot ##

ax = fig.add_subplot(2, 2, 4, projection='3d')
plt.title(r'$\omega(x,y)$')
data = np.genfromtxt('../plot_data/plot_plot_omega')
x = data[:,0]
y = data[:,1]
z = data[:,2]

xi = np.linspace(min(x), max(x))
yi = np.linspace(min(y), max(y))

X, Y = np.meshgrid(xi, yi)
Z = griddata(x, y, z, xi, yi)

surf = ax.plot_surface(X, Y, Z, rstride=2, cstride=2, cmap=cm.jet,
                       linewidth=0.02, antialiased=True)

ax.set_zlim3d(np.min(Z), np.max(Z))
fig.colorbar(surf)

plt.show()
