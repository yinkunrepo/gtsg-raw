#!/usr/bin/env python
# extract data, convert unit, and split file

import numpy as np

# input data
# ---------------------------------
input_file = 'input.txt'
volume_basis = 9.800 # A^3/atom
volume_exp = 47.014 / 5 # A^3/atom
pressure_emp = 0 # GPa, pressure correction
# ---------------------------------
# read data
dat = np.genfromtxt(input_file)

# Get middle temperature
max_T = max(dat[:,4])
min_T = min(dat[:,4])
mid_T = (max_T + min_T)*0.5*1e-3
f = open('mid_T.txt','w')
print("%.3f"%(mid_T), file=f)
f.close()

# Get unique a/a0 lattice ratios
ratio_column = dat[:,0]
lattice_ratio = np.unique(ratio_column)[::-1]

file_v = 'V.txt'
f0 = open(file_v, 'w')
print("# V(*1e-29 m^3/atom)", file=f0)
nline = dat.shape[0]
j = 0
for r in lattice_ratio:
    j = j + 1
    file_tup = 'V' + str(j) + '.txt'
    f = open(file_tup, 'w')
    print("# T(*1e+3 K), U(*1e-19 J/atom), p(*1e+10 Pa)", file=f)
    this_volume = volume_basis*(r**3)*1.0e-1 # A^3/atom ==> m^3/atom*1e-29
    print("%f" %(this_volume), file=f0)
    T = []
    p = []
    U = []
    for i in range(nline):
        this_ratio = dat[i,0]
        this_pressure = dat[i,2]*1.0e-2 # kBar ==> Pa*1e+10
        pressure_correction = pressure_emp*1.0e-1 # GPa ==> Pa*1e+10
        this_pressure = this_pressure + pressure_correction
        this_energy = dat[i,3]*1.602176634e-19*1.0e+19 # eV/atom ==> J/atom*1e-19
        delta_volume = this_volume - volume_exp*1.0e-1 # m^3/atom*1e-29
        energy_correction = -pressure_correction*delta_volume # J/atom*1e-19
        this_energy = this_energy + energy_correction
        this_temperature = dat[i,4]*1.0e-3 # K ==> K*1e+3
        if (r == this_ratio):
            T.append(this_temperature)
            p.append(this_pressure)
            U.append(this_energy)
            print("%f %f %f" %(this_temperature,this_energy,this_pressure), file=f)
    f.close()
f0.close()
