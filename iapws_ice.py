from sympy import log

from algebraic_manipulations import *

from iapws_boundaries import p_t, T_t

def gibbs_ice_I(T,p):
    t1 = 3.68017112855051e-2+5.10878114959572e-2j
    t2 = 3.37315741065416e-1+3.35449415919309e-1j
    r1 = 4.47050716285388e1+  6.56876847463481e1j
    r2 = [ -7.25974574329220e01   + -7.81008427112870e1j,
           -5.57107698030123e-5   + 4.64578634580806e-5j,
            2.34801409215913e-11+ -2.85651142904972e-11j]
    p_i0 = 1.01325e5 / p_t
    s0 = -3.32733756492168e3
    #s0 = 1.8913e2 #for absolute

    G0 = [ -6.32020233335886e5,  6.55022213658955e-1,
           -1.89369929326131e-8,  3.39746123271053e-15,
           -5.56464869058991e-22 ]
    p_i = p/p_t
    T_i = T/T_t
    
    g0_of_p = sum([ gk *(p_i-p_i0)**k for k,gk  in enumerate(G0) ])
    r2_of_p = sum([ r2k*(p_i-p_i0)**k for k,r2k in enumerate(r2) ])
    
    term_r1 = r1     *( (t1-T_i)*log(t1-T_i) + (t1+T_i)*log(t1+T_i)-2*t1*log(t1)-T_i**2/t1 )
    term_r2 = r2_of_p*( (t2-T_i)*log(t2-T_i) + (t2+T_i)*log(t2+T_i)-2*t2*log(t2)-T_i**2/t2 )
    return (g0_of_p - s0 * T + T_t*(term_r1 + term_r2))

density_ice_I,enthalpy_ice_I = density_enthalpy(gibbs_ice_I)
