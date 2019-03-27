# Code used to check numericaly the expression of theta_n(z)

import numpy as np
from zeta_toolkite import *

x, y = 0.4, -2
z = x + 1j*y
n = 14

zeta  = Zeta(z)
zeta_2 = Zeta(2*x)


def sum_cos(zeta, n):

    ret = 0

    for i in range(1, n+2):

        ret += np.cos(np.abs(zeta.z.imag)*np.log((i+1)/i) + zeta.get_gamma(i))

    return ret

#print(sum_cos(zeta, n))
#print(zeta_2.get_zeta(n+2))

print(zeta.get_gamma(n))
print(zeta.get_gamma(n))

sin_n2 = np.sin(-zeta.z.imag*np.log((n+2)/(n+1)) + zeta.get_gamma(n))

cos_n2 = np.cos(- zeta.z.imag*np.log((n+2)/(n+1)) + zeta.get_gamma(n))

sin_n1 = np.sin(-zeta.z.imag*np.log((n+1)/(n)) + zeta.get_gamma(n-1))

frac = sin_n1/(np.sin(zeta.get_theta(n))*(n+1)**x)
print(frac)
print(np.abs(zeta.get_zeta(n+1)))

sqrt = np.sqrt(frac**2 + (n+2)**(-2*zeta.z.real) + 2*cos_n2*frac*(n+2)**(-zeta.z.real))
#print(np.abs(zeta.get_zeta(n+2)))

sin_theta_n = sin_n2/(sqrt*(n+1)**x)

#print(sin_theta_n1)
#print(np.sin(zeta.get_theta(n+1)))
#print(np.sin(zeta.get_theta(n)))
