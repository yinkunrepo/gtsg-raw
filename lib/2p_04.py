#!/usr/bin/env python
# find intersection of P-T curves at constant G

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# Specify isovalue levels of Gibbs free energy
n_level = 13
G_min = -8.2
G_max = 2.2
level = np.linspace(G_min,G_max,n_level)

x_label = r'Temperature [$\times 10^{3}$ K]'
y_label = r'Pressure [$\times 10^{10}$ Pa]'

def get_cross_point(level_str):
    dat_solid = np.genfromtxt('solid_level' + level_str + '.txt')
    dat_liquid = np.genfromtxt('liquid_level' + level_str + '.txt')
    order = 3
    # solid
    x = dat_solid[:,0]
    y = dat_solid[:,1]
    pa = np.polyfit(x,y,order)
    plt.clf()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.plot(x,y,'r+')
    min_solid = min(x)
    max_solid = max(x)
    # liquid
    x = dat_liquid[:,0]
    y = dat_liquid[:,1]
    pb = np.polyfit(x,y,order)
    plt.plot(x,y,'gx')
    min_liquid = min(x)
    max_liquid = max(x)
    # look for intersection
    lower_bound = max(min_solid, min_liquid)
    upper_bound = min(max_solid, max_liquid)
    xs = np.linspace(lower_bound, upper_bound, 50)
    ys = np.polyval(pa, xs)
    plt.plot(xs, ys, 'r-')
    ys = np.polyval(pb, xs)
    plt.plot(xs, ys, 'g-')
    guess = 0.5*(lower_bound + upper_bound)
    f = np.poly1d(pa-pb)
    root = fsolve(f, guess)
    temperature = root[0]
    pressure_a = np.polyval(pa, temperature)
    pressure_b = np.polyval(pb, temperature)
    pressure = 0.5*(pressure_a + pressure_b)
    print("%s %.3f %.3f" %(level_str, pressure, temperature))
    plt.plot(temperature,pressure,'bo')    
    plt.savefig('gibbs_energy'+level_str+'.pdf',format='pdf')

print("# G(*1e-19 J/atom), P(*1e+10 Pa), T(*1e+3 K)")
for level_value in level:
    level_value_str = "%+.3f" %(level_value)
    get_cross_point(level_value_str)
