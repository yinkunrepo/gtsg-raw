#!/usr/bin/env python
# entropy calibration

import numpy as np
from gibbs_surf import *
from scipy import constants

N_A = constants.N_A

# calibrated state point
T_ref = 5.425
p_ref = 4.80

# solid phase
dat_abc = np.genfromtxt('solid_dat_abc.txt')
dat_eos = np.genfromtxt('solid_dat_eos.txt')
[V,S,p,T,U,F,H,G]=potential_G(T_ref,p_ref,dat_abc,dat_eos)
G_solid = G

# liquid phase
dat_abc = np.genfromtxt('liquid_dat_abc.txt')
dat_eos = np.genfromtxt('liquid_dat_eos.txt')
[V,S,p,T,U,F,H,G]=potential_G(T_ref,p_ref,dat_abc,dat_eos)
V_liquid = V
S_liquid = S
p_liquid = p
T_liquid = T
U_liquid = U

S_before_shift = S_liquid
S_after_shift = (U_liquid + p_liquid*V_liquid - G_solid) / T_liquid
delta_sx = S_after_shift - S_before_shift
f = open('liquid_delta_S.txt','w')
print("%.6f" %(delta_sx), file=f)
f.close()
