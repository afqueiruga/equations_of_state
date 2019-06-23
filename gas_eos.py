import numpy as np
import sympy
from sympy import log, exp

"""
Fits for gases
"""
from iapws_boundaries import R

properties = {'CH4':{
    'molar_weight':1.6043e-2,
     'T_crit':1.9056e2,
     'p_crit':4.600155e6,
     'V_crit':9.9e-2,
     'T_crit':2.880e-1, 
     'omega':1.10002,
     'A':[4.568,-8.9750e-3,  3.631e-5, -3.407e-8,  1.091e-11] 
     }
}
@np.vectorize
def pressure_peng_robinson(rho,T, gas='CH4'):
    "Peng Robinson is closed form for pressure"
    T_c   = properties[gas]['T_crit'] # critical temperature
    p_c   = properties[gas]['p_crit'] # critical pressure
    omega = properties[gas]['omega']  # acentric factor
    Mw = properties[gas]['molar_weight']
    Vm = Mw / rho
    a = 0.45724*R**2*T_c**2/p_c
    b = 0.07780*R*T_c/p_c
    kappa = 0.37464+1.54226*omega-0.26992*omega**2
    alpha = (1.0+kappa*(1.0-(T/T_c)**0.5))**2

    p = R*T/(Vm-b) - a*alpha/(Vm**2+2*b-b**2) 
    return p


def rho_peng_robinson(p, T, gas='CH4'):
    "Invert peng robinsion to get rho"
    rho = sympy.Symbol('rho')
    expr = pressure_peng_robinson(rho,T)-p
    return float(sympy.nsolve(expr,rho,10.0))

def enthalpy_gas(p,T):
    pass
