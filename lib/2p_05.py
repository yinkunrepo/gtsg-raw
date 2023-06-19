#!/usr/bin/env python
# fit the Kechin equation

import numpy as np
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

# pressure range and step (in 1e+10 GPa)
P_min = 1.3
P_max = 13.6
P_step = 0.1
# let T0 be fixed (in 1e+3 K)
T0 = 4.020
###

def kechin_equation(x, P0, a, b, c):
    return T0*(1.0 + (x-P0)/a)**b*np.exp(-(x-P0)/c)

dat = np.genfromtxt('cross_points.txt')
pressures = dat[:,1]
temperatures = dat[:,2]

guess = np.array([4.42, 4.88, 0.413, 86.0])
popt, pcov = curve_fit(kechin_equation, pressures, temperatures, p0=guess, maxfev=1000000)

# write down the parameters of Kechin equation with units
f = open('kechin_fit_params_with_units.txt','w')
print("T0 = %.0f (K)" %(T0*1000), file=f)
print("P0 = %.1f (GPa)" %(popt[0]*10), file=f)
print("a = %.1f (GPa)" %(popt[1]*10), file=f)
print("b = %.3f (dimensionless)" %(popt[2]), file=f)
print("c = %.2e (GPa)" %(popt[3]*10), file=f)
f.close()

# write without units in a single line
f = open('kechin_fit_params.txt','w')
print("%.3f %.2f %.2f %.3f %.2e" %(T0, popt[0], popt[1], popt[2], popt[3]), file=f)
f.close()

f = open('fit_kechin_data.txt','w')
xs = np.arange(P_min, P_max+P_step, P_step)
ns = len(xs)
ys = np.empty(ns)
params = np.genfromtxt('kechin_fit_params.txt')
print("# P(1e+10Pa), T(1e+3K)", file=f)
for i in range(ns):
    ys[i] = kechin_equation(xs[i], params[1], params[2], params[3], params[4])
    print("%.1f %.3f" %(xs[i],ys[i]), file=f)
f.close()

import matplotlib.pyplot as plt
x_label = r'Temperature [$\times 10^{3}$ K]'
y_label = r'Pressure [$\times 10^{10}$ Pa]'
plt.xlabel(y_label)
plt.ylabel(x_label)
plt.xlim(1.3,13.6)
plt.plot(pressures,temperatures,'ko',mfc='none')
plt.plot(xs, ys, 'r-')
plt.savefig('kechin_fit.pdf',format='pdf')
