#!/usr/bin/env python
# fit linear functions U(T)=a*T+b and p(T)=d*T+e

import numpy as np

dat = np.genfromtxt('V.txt')
file_out = 'abde.txt'
f = open(file_out,'w')
n = len(dat)
for i in range(n):
    file_tup = 'V' + str(i+1) + '.txt'
    dat = np.genfromtxt(file_tup)
    T = dat[:,0]
    U = dat[:,1]
    p = dat[:,2]
    [a,b] = np.polyfit(T,U,1)
    [d,e] = np.polyfit(T,p,1)
    print("%f %f %f %f" %(a,b,d,e), file=f)
f.close()
