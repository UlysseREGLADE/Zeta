import numpy as np
from mpmath import zeta, gamma, sin
import matplotlib.pyplot as plt

N_max = 300
z = 0.5+14.1347251417346937904572519835624766j

def b(N, t):

    if(N==t):
        return 1

    # mpmath_sin = sin((t-N)*np.pi)
    # float_sin = float(mpmath_sin.real)+float(mpmath_sin.imag)*1j
    #
    # return float_sin/((t-N)*np.pi)

    try:
        mpmath_gamma = gamma(-t)
        float_gamma = float(mpmath_gamma.real)+float(mpmath_gamma.imag)*1j
        ret = (((-1)**N)/((N-t)*gamma(N+1)))/float_gamma
    except:
        ret = 0

    return ret

T = np.linspace(0, 10, 201)
res = np.zeros(len(T), dtype=complex)
zeta_0_n = 0

plt.figure()

for N in range(N_max):

    for i in range(len(T)):
        t = T[i]
        res[i] += b(N, t)*zeta_0_n

    zeta_0_n += (N+1)**(-z)

    if(N <= np.max(T)):
        plt.plot(zeta_0_n.real, zeta_0_n.imag, 'bo')

plt.plot(res.real, res.imag)
#plt.plot(T, T*(T+1)/2)
plt.show()
