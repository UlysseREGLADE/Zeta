import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from mpmath import zeta

z = 0.5 + 14.1347251417346937904572519835624766j
#z = 0.5 + 21.0220396387715549926284795938969162j
#z = 0.5 + 25.0108575801456887632137909925627734j
#z = 0.5 + 30.4248761258595132103118975305839571j
#z = 0.5 + 32.935061587739189690662368964074741j
#z = 0.9 + 10j
N = 10000

T = np.linspace(-0.85, 0, 850+1)
T_r = np.array([-1+1/6, -1+1/4, -1+1/3 , -0.5, -1/3, -1/4, -1/6])
Q, R = (T//1).astype(int), T%1
Q_r, R_r = (T_r//1).astype(int), T_r%1
partial_sum = np.zeros((len(T)), dtype=complex)
partial_sum_r = np.zeros((len(T_r)), dtype=complex)

zeta_z = 0
for n in range(1, N):
    zeta_z += n**-z

def cont_partial_sum(q, r):

    ret = 0
    for n in range(q+1, N):
        ret += (n+r)**-z

    return zeta_z - ret

for i in range(len(T)):

    partial_sum[i] = cont_partial_sum(Q[i], R[i])

for i in range(len(T_r)):

    partial_sum_r[i] = cont_partial_sum(Q_r[i], R_r[i])

print(np.abs(partial_sum_r[0] - partial_sum_r[-1]))
print(np.abs(partial_sum_r[1] - partial_sum_r[-2]))
print(np.abs(partial_sum_r[2] - partial_sum_r[-3]))

plt.figure()
#plt.plot(T, partial_sum)
plt.plot(partial_sum.real, partial_sum.imag)
plt.plot(partial_sum_r.real, partial_sum_r.imag, "o")
plt.plot(zeta(z).real, zeta(z).imag, "or")
#plt.plot([0, T[-1]], [zeta_z, zeta_z])
plt.axis('equal')
plt.show()
