#!/usr/bin/env python
# specify the properties V, T, S at a reference state,
# derive the preliminary coefficients of U(S) at constant V,
# generate coefficients of U(V) at constant S, i.e., dat_eos.txt,
# generate coefficients of U(S) at constant V again, i.e., dat_abc.txt,
# to make the surface more smooth.

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import leastsq

# whether output full information and intermediate graphics to debug
debug = False

# load data
V = np.genfromtxt('V.txt')
nV = len(V)
dat = np.genfromtxt('abde.txt')
a = dat[:,0]
b = dat[:,1]
d = dat[:,2]
e = dat[:,3]

# set properties at the reference state
refS = np.genfromtxt('refS.txt')
ref_index = int(np.genfromtxt('ref_index.txt'))
refV = V[ref_index]
f = open('refV.txt','w')
print("%f"%(refV), file=f)
f.close()
refT = np.genfromtxt('refT.txt')
refU = a[ref_index]*refT + b[ref_index]

# get V-p data at the reference T
p = np.empty(nV)
for i in range(nV):
    p[i] = d[i]*refT + e[i]

# fit V-p to Equation of State (EoS)
# here we use a 3rd-order Birch-Murnaghan EoS
def eos(V,B0,BP,V0):
    eta = (V0/V)**(1./3.)
    p = 1.5*B0*(eta**7-eta**5)*(1+3./4.*(BP-4)*(eta**2-1))
    return p

popt, pcov = curve_fit(eos, V,p, maxfev=1000000)
B0 = popt[0]
BP = popt[1]
V0 = popt[2]

# derive preliminary coefficient c
# get Helmholtz free energy by integration of EoS (without constant A0)
def A_without_constant(V,B0,BP,V0):
    eta = (V0/V)**(1./3.)
    return 9.*V0*B0/16.*((eta**2-1)**3*BP+(eta**2-1)**2*(6-4*eta**2))

file_out = 'preliminary_c.txt'
f = open(file_out,'w')
for i in range(nV):
    thisV = V[i]
    thisU = a[i]*refT + b[i]
    thisA = A_without_constant(thisV,B0,BP,V0)
    refA = A_without_constant(refV,B0,BP,V0)
    thisS = ((thisU-refU+refT*refS) - (thisA-refA))/refT
    c = a[i]*refT/np.exp(thisS/a[i])
    print("%f" %(c), file=f)
f.close()

# generate coefficients of U(V) at constant S, i.e., dat_eos.txt
c = np.genfromtxt('preliminary_c.txt')
def func_ST(T,a,c):
    S = a*np.log(a*T/c)
    return S

def func_TS(S, a, c):
    T = c/a*np.exp(S/a)
    return T

# determine the entropy range to cover
S_all = []
for i in range(nV):
    file_tup = 'v' + str(i+1) + '.txt'
    dat = np.genfromtxt(file_tup)
    T = dat[:,0]
    for thisT in T:
        thisS = func_ST(thisT,a[i],c[i])
        S_all.append(thisS)

minS = min(S_all)
maxS = max(S_all)

if debug:
    print("entropy range:")
    print("minimum = %f, maximum = %f"%(minS,maxS))
    for i in range(nV):
        mintemp = func_TS(minS,a[i],c[i])
        maxtemp = func_TS(maxS,a[i],c[i])
        print("V = %f, minimum T = %f, maximum T = %f"%(V[i],mintemp,maxtemp))
        
def func_UV(V,E0,B0,BP,V0):
    eta = (V0/V)**(1./3.)
    U = E0 + 9.*B0*V0/16.*((eta**2-1)**3*BP + (eta**2-1)**2*(6-4*eta**2))
    return U

def bm3(parameters, V):
    E0, B0, BP, V0 = parameters
    eta = (V0/V)**(1./3.)
    U = E0 + 9.*B0*V0/16.*((eta**2-1)**3*BP + (eta**2-1)**2*(6-4*eta**2))
    return U

def objective(pars, y, x):
    err = y - bm3(pars, x)
    return err

def func_US(S,a,b,c):
    U = c*np.exp(S/a)+b
    return U

file_out = 'dat_eos.txt'
f = open(file_out, 'w')
nS = 10 # number of entropies to sample
S = np.linspace(minS, maxS, nS)

def parabola(x,a,b,c):
   return a + b*x + c*x**2

fig = plt.figure()
plt.subplot(111)
for i in range(nS):
    thisS = S[i]
    U = np.empty(nV)
    for j in range(nV):
        U[j] = func_US(thisS,a[j],b[j],c[j])
    # plot figure
    plt.xlabel('V (1.0e-29 m^3/atom)')
    plt.ylabel('U (1.0e-19 J/atom)')
    plt.title('Constant S = ' + "%.3f" %(thisS) + '*1.0e-22 J/K/atom')
    plt.plot(V,U,'o')
    # estimate the initial guess
    p0 = [min(U),1,1]
    popt, pcov = curve_fit(parabola,V,U,p0,maxfev=1000000)
    parabola_parameters = popt
    a0 = parabola_parameters[0]
    b0 = parabola_parameters[1]
    c0 = parabola_parameters[2]
    parabola_vmin = -b0/2/c0
    E00 = parabola(parabola_vmin,a0,b0,c0)
    B00 = 2*c0*parabola_vmin
    BP0 = 4
    initial_guess = [E00,B00,BP0,parabola_vmin]
    if i == 0:
       guess = initial_guess
    plsq = leastsq(objective, guess, args=(U,V), maxfev=1000000)
    popt = plsq[0]
    E0 = popt[0]
    B0 = popt[1]
    BP = popt[2]
    V0 = popt[3]
    guess = [E0,B0,BP,V0]
    print("%f %f %f %f %f" %(thisS, E0,B0,BP,V0), file=f)
    if debug:
        print("%f %f %f %f %f" %(thisS, E0,B0,BP,V0))
    nk = 50
    V_s = np.linspace(min(V),max(V),nk)
    U_s = np.empty(nk)
    for k in range(nk):
        U_s[k] = func_UV(V_s[k],E0,B0,BP,V0)
    if debug:
        plt.plot(V_s,U_s,'k-')
        plt.show()
f.close()

# generate coefficients of U(S) at constant V again, i.e., dat_abc.txt
# use the first line of preliminary a,b,c as the guess
dat_abde = np.genfromtxt('abde.txt')
dat_c = np.genfromtxt('preliminary_c.txt')
guess_a = dat_abde[0,0]
guess_b = dat_abde[0,1]
guess_c = dat_c[0]
guess = [guess_a, guess_b, guess_c]

file_out = 'dat_abc.txt'
f = open(file_out, 'w')
dat = np.genfromtxt('dat_eos.txt')
S = dat[:,0]
nS = dat.shape[0]
for i in range(nV):
    thisV = V[i]
    U = np.empty(nS)
    for j in range(nS):
        E0 = dat[j,1]
        B0 = dat[j,2]
        BP = dat[j,3]
        V0 = dat[j,4]
        U[j] = func_UV(thisV,E0,B0,BP,V0)
    popt, pcov = curve_fit(func_US,S,U,guess,maxfev=1000000)
    a = popt[0]
    b = popt[1]
    c = popt[2]
    print("%f %f %f %f" %(thisV,a,b,c), file=f)
f.close()
