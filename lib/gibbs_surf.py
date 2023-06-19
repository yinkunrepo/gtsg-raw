import numpy as np
from scipy.optimize import curve_fit, fsolve
import warnings
warnings.filterwarnings('ignore')

def func_US(S, a, b, c):
    U = c*np.exp(S/a) + b
    return U

def func_TS(S, a, c):
    T = c/a*np.exp(S/a)
    return T

def func_ST(T, a, c):
    S = a*np.log(a*T/c)
    return S

# BirchMurnaghan equation from PRB 70, 224107
def func_UV(V, E0, B0, BP, V0):
    eta = (V0/V)**(1./3.)
    U = E0 + 9.*B0*V0/16.*( (eta**2-1)**3*BP + (eta**2-1)**2*(6-4*eta**2) )
    return U

def func_pV(V, B0, BP, V0):
    eta = (V0/V)**(1./3.)
    p = 1.5*B0*(eta**7-eta**5)*(1+3./4.*(BP-4)*(eta**2-1))
    return p

def func_Vp(p, B0, BP, V0, dat_abc):
    minV = min(dat_abc[:,0])
    maxV = max(dat_abc[:,0])
    def f(x):
        eta = (V0/x)**(1./3.)
        return 1.5*B0*(eta**7-eta**5)*(1+3./4.*(BP-4)*(eta**2-1)) - p
    # determine the guess volume
    n_V = 50
    V_list = np.linspace(minV, maxV, n_V)
    delta_p_list = np.empty(n_V)
    for j in range(n_V):
        volume = V_list[j]
        pressure = func_pV(volume, B0, BP, V0)
        delta_p_list[j] = abs(pressure - p)
    min_delta_idx = np.argmin(delta_p_list)
    V_guess = V_list[min_delta_idx]
    result = fsolve(f, V_guess)
    V = result[0]
    return V
# end of BirchMurnaghan equation

def fit_US(entropies, energies, dat_abc):
    guess_abc = dat_abc[0,1:]
    popt, pcov = curve_fit(func_US, entropies, energies, guess_abc, maxfev=1000000)
    a = popt[0]
    b = popt[1]
    c = popt[2]
    return [a, b, c]

def fit_UV(volumes, energies, dat_eos):
    guess_eos = dat_eos[0,1:]
    popt, pcov = curve_fit(func_UV, volumes, energies, guess_eos, maxfev=1000000)
    E0 = popt[0]
    B0 = popt[1]
    BP = popt[2]
    V0 = popt[3]
    return [E0, B0, BP, V0]

def at_const_S(S, dat_abc, dat_eos):
    dat = dat_abc
    nd = dat.shape[0]
    volumes = dat[:,0]
    energies = np.empty(nd)
    for i in range(nd):
        a = dat[i,1]
        b = dat[i,2]
        c = dat[i,3]
        energies[i] = func_US(S,a,b,c)
    [E0, B0, BP, V0] = fit_UV(volumes, energies, dat_eos)
    return [E0,B0,BP,V0]

def at_const_V(V, dat_abc, dat_eos):
    dat = dat_eos
    nd = dat.shape[0]
    entropies = dat[:,0]
    energies = np.empty(nd)
    for i in range(nd):
        E0 = dat[i,1]
        B0 = dat[i,2]
        BP = dat[i,3]
        V0 = dat[i,4]
        energies[i] = func_UV(V, E0,B0,BP,V0)
    [a,b,c] = fit_US(entropies,energies,dat_abc)
    return [a,b,c]

def potential_U(S, V, dat_abc, dat_eos):
    [E0,B0,BP,V0] = at_const_S(S,dat_abc,dat_eos)
    p = func_pV(V,B0,BP,V0)
    U1 = func_UV(V,E0,B0,BP,V0)
    [a, b, c] = at_const_V(V,dat_abc,dat_eos)
    T = func_TS(S, a, c)
    U2 = func_US(S, a, b, c)
    U = (U1 + U2)*0.5
    A = U - T*S
    H = U + p*V
    G = U - T*S + p*V
    return [V, S, p, T, U, A, H, G]

def potential_A(T, V, dat_abc, dat_eos):
    [a, b, c] = at_const_V(V,dat_abc,dat_eos)
    S = func_ST(T, a, c)
    [V, S, p, T, U, A, H, G] = potential_U(S, V, dat_abc,dat_eos)
    return [V, S, p, T, U, A, H, G]

def potential_H(S, p, dat_abc, dat_eos):
    [E0,B0,BP,V0] = at_const_S(S,dat_abc,dat_eos)
    V = func_Vp(p,B0,BP,V0,dat_abc)
    [V, S, p, T, U, A, H, G] = potential_U(S, V,dat_abc,dat_eos)
    return [V, S, p, T, U, A, H, G]

def potential_G(T, p, dat_abc, dat_eos):
    temperature = T
    pressure = p
    nd = dat_eos.shape[0]
    entropies = dat_eos[:,0]
    enthalpies = np.empty(nd)
    for i in range(nd):
        [V,S,p,T,U,A,H,G] = potential_H(entropies[i],pressure,dat_abc,dat_eos)
        enthalpies[i] = H
    [a,b,c] = fit_US(entropies,enthalpies,dat_abc)
    entropy = func_ST(temperature,a,c)
    [V,S,p,T,U,A,H,G] = potential_H(entropy,pressure,dat_abc,dat_eos)
    return [V, S, p, T, U, A, H, G]
