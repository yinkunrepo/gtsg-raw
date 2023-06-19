#!/usr/bin/env python
# get isovalue lines of Gibbs free energy from contour plot
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

# Specify isovalue levels of Gibbs free energy
n_level = 13
G_min = -8.2
G_max = 2.2
level = np.linspace(G_min,G_max,n_level)

# input data
dat_sol = np.genfromtxt('data_solid.txt')
dat_liq = np.genfromtxt('data_liquid.txt')

# Solid phase
# X=>temperature, Y=>pressure, Z=>Gibbs free energy
X = dat_sol[:,3]
Y = dat_sol[:,2]
Z = dat_sol[:,7]
xi = np.linspace(X.min(),X.max(),60)
yi = np.linspace(Y.min(),Y.max(),60)
zi = griddata((X,Y), Z, (xi[None,:], yi[:,None]), method='cubic')
xig, yig = np.meshgrid(xi, yi)
cset = plt.contour(xig, yig, zi, level)
for i in range(n_level):
    level_value = cset.levels[i]
    level_value_str = "%+.3f" %(level_value)
    file_name = "solid_level" + level_value_str + ".txt"
    path = cset.collections[i].get_paths()[0]
    vert = path.vertices
    f = open(file_name, 'w')
    for j in range(vert.shape[0]):
        print("%f %f" %(vert[j,0],vert[j,1]), file=f)
    f.close()

# Liquid phase
# X=>temperature, Y=>pressure, Z=>Gibbs free energy
X = dat_liq[:,3]
Y = dat_liq[:,2]
Z = dat_liq[:,7]
xi = np.linspace(X.min(),X.max(),60)
yi = np.linspace(Y.min(),Y.max(),60)
zi = griddata((X,Y), Z, (xi[None,:], yi[:,None]), method='cubic')
xig, yig = np.meshgrid(xi, yi)
cset = plt.contour(xig, yig, zi, level)
for i in range(n_level):
    level_value = cset.levels[i]
    level_value_str = "%+.3f" %(level_value)
    file_name = "liquid_level" + level_value_str + ".txt"
    path = cset.collections[i].get_paths()[0]
    vert = path.vertices
    f = open(file_name, 'w')
    for j in range(vert.shape[0]):
        print("%f %f" %(vert[j,0],vert[j,1]), file=f)
    f.close()
