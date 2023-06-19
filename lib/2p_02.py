#!/usr/bin/env python
# generate grid data
import numpy as np
from gibbs_surf import *

# Specify the mesh grid of data
nV = 40
nT = 40

# Specify the maximum and minimum values of temperature to cover
T_min = 2.3
T_max = 8.0

def generate_grid(phase='solid'):
    dat_abc = np.genfromtxt(phase+'_dat_abc.txt')
    dat_eos = np.genfromtxt(phase+'_dat_eos.txt')
    volumes = dat_abc[:,0]
    V_min = min(volumes)
    V_max = max(volumes)
    f = open('data_'+phase+'.txt', 'w')
    V_data = np.linspace(V_min, V_max, nV)
    T_data = np.linspace(T_min, T_max, nT)
    for i in range(nV):
        for j in range(nT):
            [V,S,P,T,E,A,H,G] = potential_A(T_data[j],V_data[i],dat_abc,dat_eos)
            print("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f" %(V,S,P,T,E,A,H,G), file=f)
    f.close()

generate_grid(phase='solid')
generate_grid(phase='liquid')
