import sympy as sy
from sympy.utilities.lambdify import lambdify
import mpmath as mp
import matplotlib.pyplot as plt
import numpy as np

M_max = 6

M_max += 1

N, z = sy.symbols('N z')
x = sy.symbols('x')

R = [sy.S('0')]

for M in range(M_max):

    new_R = R[-1]
    new_R += (R[-1].subs(N, N+1)-R[-1]+(N+1)**(-z))*N/(M+z-1)

    for i in range(M_max):

        I = i+1

        new_R = new_R.subs((N+I)**-z, x**I)

    new_R = sy.apart(new_R, x)

    for i in range(M_max):

        I = M_max - i

        new_R = new_R.subs(x**I, (N+I)**-z)

    R.append(new_R)

    print(M)

f = []
for M in range(M_max):
    f.append(lambdify((N, z), R[M]))

N = np.linspace(0, M_max, 50*(M_max)+1)
z=2

#f_N = -f(N, z)
#f_N += float(mp.zeta(z).real)+float(mp.zeta(z).imag)*1j

#N = np.linspace(-0.75, 0, 101)
plt.figure()

for M in range(M_max):

    f_N_neg = -f[M](N, z)
    f_N_neg = f_N_neg + np.array(mp.zeta(z).real, dtype=float)
    print(N.shape)
    print(f_N_neg.shape)
    plt.plot(N, f_N_neg.real, label = M)

#plt.plot(f_N.real, f_N.imag)
zeta = 1
for n in range(2, M_max-1):
    plt.scatter(zeta.real, zeta.imag, color='r')
    zeta += n**-z
plt.scatter(mp.zeta(z).real,mp.zeta(z).imag, color='b')
plt.legend()
plt.show()
