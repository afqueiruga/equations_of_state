import numpy as np
from numpy import exp

"""
Boundaries between the phases of water, taken from the IAPWS 1995 release.

http://www.iapws.org/relguide/IAPWS-95.html
"""

R = 0.461526e3

triple_point   = T_t,p_t = (273.16, 611.657)
critical_point = T_c,p_c = (647.096, 22.064e6)
rho_critical   = rho_c = 322.0

def vapor_pressure(T):
    a = [ None, # Skip to follow IAPWS indexing
        -7.85951783,
         1.84408259,
        -11.7866497,
         22.6807411,
        -15.9618719,
         1.80122502
    ]
    th = 1.0 - T/T_c
    ex = T_c/T * (a[1]*th + a[2]*th**1.5 + a[3]*th**3.0
                  + a[4]*th**3.5 + a[5]*th**4 + a[6]*th**7.5 )
    return p_c * exp(ex)

melting_pressure_range_I = (251.165,273.16)
def melting_pressure_I(T):
    T_n = melting_pressure_range_I[1]
    p_n = 0.000611657e6
    th = T/T_n
    return p_n*(1.0-0.626e6*(1.0-th**-3)+0.197135e6*(1.0-th**21.2))
melting_pressure_range_III = (251.165,256.164)
def melting_pressure_III(T):
    T_n = melting_pressure_range_III[0]
    p_n = 209.9e6
    th = T/T_n
    return p_n*( 1.0-0.295252*(1.0-th**60) )
melting_pressure_range_V = (256.164,273.31)
def melting_pressure_V(T):
    T_n = melting_pressure_range_V[0]
    p_n = 350.1e6
    th = T/T_n
    return p_n*( 1.0 - 1.18721*(1.0-th**8) )

melting_pressure_range_VI = (273.31, 355.0)
def melting_pressure_VI(T):
    T_n = melting_pressure_range_VI[0]
    p_n = 632.4e6
    th = T/T_n
    return p_n * ( 1.0 - 1.07476*(1.0-th**4.6) )
melting_pressure_range_VII = (355.0, 715.0)
def melting_pressure_VII(T):
    T_n = melting_pressure_range_VII[0]
    p_n = 2216.0e6
    th = T/T_n
    return p_n * exp( 1.73683*(1.0-th**-1) - 0.0544606*(1.0-th**5)
                     + 0.806106e-7*(1.0-th**22) )

sublimation_pressure_range = (100, T_t)
def sublimation_pressure(T):
    p_n = 0.000611657e6
    th = T/T_t
    return p_n * exp( -13.928169*(1.0-th**-1.5) + 34.7078238*(1.0-th**-1.25) )

# Numerically import one
import sympy
T = sympy.Symbol('T')
P_expr = melting_pressure_I(T)
@np.vectorize
def melting_temperature_I(P):
    return float(sympy.nsolve(P_expr-P,T, 0.5*(melting_pressure_range_I[1]-melting_pressure_range_I[0]) ))


def plot_boundaries():
    """Plot the phase boundaries"""
    from matplotlib import pylab as plt
    plt.semilogy(triple_point[0],triple_point[1],'o')
    plt.semilogy(critical_point[0],critical_point[1],'o')

    for f,r in [
            (vapor_pressure,(T_t,T_c)),
            (melting_pressure_I,melting_pressure_range_I),
            (melting_pressure_III,melting_pressure_range_III),
            (melting_pressure_V,melting_pressure_range_V),
            (melting_pressure_VI,melting_pressure_range_VI),
            (melting_pressure_VII,melting_pressure_range_VII),
            (sublimation_pressure,sublimation_pressure_range)
            ]:
        _ts = np.linspace(r[0],r[1])
        plt.semilogy(_ts,f(_ts),'-')
    plt.xlim(100,900)
    plt.ylim(1,10.0e10)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Pressure (Pa)')
    plt.text(120,1e3,'Ice I')
    plt.text(200,1e10,'Ices III, V, VI, VII')
    plt.text(400,1e3,'Gas')
    plt.text(400,1e7,'Liquid')
    plt.text(700,1e8,'Supercritical')
    
