#Code used to visualize the convergence of the Riemann series

import numpy as np
from zeta_toolkite import *
import matplotlib.pyplot as plt

z = 1.5 - 3j

zeta = Zeta(z)

n = np.arange(1000)

print(zeta.get_c(0))

plt.figure()
plt.title("c_n(z) and zeta_n(z) for z = "+str(z))
plt.plot(zeta.get_zeta(n).real, zeta.get_zeta(n).imag, 'r+', label="zeta_n(z)")
plt.plot(zeta.get_c(n).real, zeta.get_c(n).imag, 'b+', label="c_n(z)")
plt.plot(zeta.get_zeta().real, zeta.get_zeta().imag, 'go', label="zeta(z)")
plt.axis('equal')
plt.legend()
plt.grid()
plt.show(True)
