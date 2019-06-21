import numpy as np
from sympy import log, exp, Piecewise

"""
Fits for methane gas and methane hydrate. 

Based on Ballard, 2002 and Moridis 2003.
See Moridis, Queiruga, and Regan, 2019.
"""
Molar_weight_ch4   = Mg = 16.04 #kg/mol
Molar_weight_water = Mw = 18.02 #kg/mol
Hydration_Number = NH = 6
Mh = NH*Mw+Mg
def enthalpy_hydrate(T,p):
    return (T-273.15)*2100.0
def density_hydrate(T,p):
    dT = T-298.15
    dP = p-1.0e5
    a1 = 3.38496e-4
    a2 = 5.40099e-7
    a3 = -4.76946e-11
    a4 = 1.0e-10
    v0 = 1000*Molar_weight_ch4 / (22.712*NH)
    return (v0*exp( a1*dT + a2*dT**2 + a3*dT**3 + a4*dP))**-1


hyd_stab_T_d = 0.0
pe_hyd_stab_A = [
    -1.94138504464560e5, 3.31018213397926e3,
    -2.25540264493806e1, 7.67559117787059e-2,
    -1.30465829788791e-4, 8.86065316687571e-8
]
pe_hyd_stab_B = [
    -4.38921173434628e01, 7.76302133739303e-1,
    -7.27291427030502e-3, 3.85413985900724e-5,
    -1.03669656828834e-7, 1.09882180475307e-10
]
pe_hyd_stab_C = [
    9.652117566301735e-1, 5.563942679876470e-2,
    2.934835672207024e-2, 7.696735279663661e-3,
    -6.147609081030884e-3, -1.931115655395969e-3,
    6.350841470341581e-4, 1.682282887657391e-4
]
def horner(x,a):
    r = a[-1]
    for i in xrange(len(a)-2,-1,-1):
        r += r*x+a[i]
    return r
@np.vectorize
def pe_hyd_stab(T):
    pa = 1.0e6*exp(horner(T,pe_hyd_stab_A))
    pb = 1.0e6*exp(horner(T,pe_hyd_stab_B))
    pc = 1.0e6*exp(horner(T-273.15,pe_hyd_stab_C))
    return Piecewise(
        (pa,T>275.0),
        (pa+2*(275.0-T)*(pc-pa), T>274.5),
        (pc, T>272.0),
        (pc+2*(272-T)*(pb-pc), T>271.5),
        (pb,True)
    )
quadruple_point = (273.2, pe_hyd_stab(273.2) )

def plot_boundaries():
    from matplotlib import pylab as plt
    _ts = np.linspace(150,320)
    plt.semilogy(_ts,pe_hyd_stab(_ts),'-')
    plt.semilogy(quadruple_point[0],quadruple_point[1],'o')
