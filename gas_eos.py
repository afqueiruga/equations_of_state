import numpy as np
from sympy import log, exp

"""
Fits for gases
"""
from iapws_boundaries import R

properties = {'CH4':{
     'T_crit':1.9056d+02,
     'p_crit':4.600155d+06,
     'omega':1.1000d-02,
     }
}
def pressure_peng_robinson(V,T, gas='CH4'):
    "Peng Robinson is closed form for pressure"
    T_c   = properties[gas]['T_crit'] # critical temperature
    p_c   = properties[gas]['p_crit'] # critical pressure
    omega = properties[gas]['omega']  # acentric factor

    a = 0.45724*R**2*T_c**2/p_c
    b = 0.07780*R*T_c/p_c
    kappa = 0.37464+1.54226*omega-0.26992*omega**2
    alpha = (1.0+kappa*(1.0-(T/T_c)**0.5))**2

    p = R*T/(V-b) - a*alpha/(V**2+2*b-b**2) 
    return p

def enthalpy_gas(p,T):
    pass
