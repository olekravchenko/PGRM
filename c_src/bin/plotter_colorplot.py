#!/usr/bin/python
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np
import sys

fig = plt.figure()
#ax = fig.gca(projection='3d')

#plt.title(sys.argv[2])

data = np.genfromtxt(sys.argv[1])
x = data[:,0]
y = data[:,1]
z = data[:,2]

xi = np.linspace(min(x), max(x))
yi = np.linspace(min(y), max(y))

X, Y = np.meshgrid(xi, yi)
Z = griddata(x, y, z, xi, yi)

plt.contourf(X, Y, Z, 100)

#ax.set_zlim3d(np.min(Z), np.max(Z))
plt.colorbar()

plt.show()
