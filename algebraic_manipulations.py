import numpy as np
import sympy

def density_enthalpy(gibbs):
    T,p = sympy.symbols('T p')
    g = gibbs(T,p)
    density = 1/g.diff(p)
    enthalpy = g - T * g.diff(T)
    rhovec = np.vectorize(sympy.lambdify([T,p],density))
    hvec = np.vectorize(sympy.lambdify([T,p],enthalpy))
                   
    return (lambda x,y : np.real(rhovec(x,y))),\
           (lambda x,y : np.real(  hvec(x,y)))

def pressure_enthalpy_from_helmholtz(f_func):
    rho,T = sympy.symbols('rho T')
    f = f_func(T,rho)
    p = rho**2 * f.diff(rho)
    h = f - T*f.diff(T) + rho*f.diff(rho)
    return np.vectorize(sympy.lambdify([T,rho],p)),\
            np.vectorize(sympy.lambdify([T,rho],h))
