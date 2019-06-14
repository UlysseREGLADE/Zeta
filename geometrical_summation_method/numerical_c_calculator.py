#Code used to compute the numerical expression of c_n(z)

import numpy as np
import sympy as sy
from zeta_toolkite import *

def A(n, z):

    alpha = sy.pi/2 - sy.atan((1-sy.re(z))/sy.im(z))

    e = sy.exp(sy.I*alpha)

    ns = sy.S(str(n))

    return sy.Matrix([[sy.im(e*(ns+1)**-z), -sy.re(e*(ns+1)**-z)],
                      [sy.im(e*(ns+2)**-z), -sy.re(e*(ns+2)**-z)]])

def zeta(n, z):

    if(n==0):
        return sy.S('0')

    zeta_n = sy.S('1')
    for i in range(2, n+1):
        zeta_n += (sy.S(str(i)))**-z

    return zeta_n

def B(n, z):

    alpha = sy.pi/2 - sy.atan((1-sy.re(z))/sy.im(z))

    ns = sy.S(str(n))

    e = sy.exp(sy.I*alpha)

    return sy.Matrix([[sy.re(zeta(ns, z))*sy.im(e*(ns+1)**-z)   -    sy.im(zeta(ns, z))*sy.re(e*(ns+1)**-z)],
                      [sy.re(zeta(ns+1, z))*sy.im(e*(ns+2)**-z) - sy.im(zeta(ns+1, z))*sy.re(e*(ns+2)**-z)]])

def c(n, z):

    A0 = sy.simplify(A(n, z))
    invA0 = sy.simplify(A0**-1)
    B0 = sy.simplify(B(n, z))

    return sy.simplify(invA0*B0)

x = sy.symbols('x', real=True)
y = sy.symbols('y', real=True, negative=True)

z = x+sy.I*y

if(__name__=="__main__"):

    n = 0

    A0 = sy.simplify(A(n, z))

    invA0 = sy.simplify(A0**-1)

    print("invA0")
    sy.pprint(invA0)
    print("\n")

    B0 = sy.simplify(B(n, z))

    print("B0")
    sy.pprint(B0)
    print("\n")

    c_n = sy.simplify(invA0*B0)

    print("c_"+str(n))
    sy.pprint(c_n)
    print("\n")


    c_n_x = c_n[0]
    c_n_y = c_n[1]

    lambda_x = sy.lambdify((x, y), c_n_x)
    lambda_y = sy.lambdify((x, y), c_n_y)

    zn = 0.5+2j

    zeta = Zeta(zn)

    lamba = np.exp(-zn*np.log(2))*np.exp(1j*np.arctan((zn.real-1)/zn.imag))
    ps = np.conj(1-zn).real*lamba.real + np.conj(1-zn).imag*lamba.imag
    c_0_simplified = np.conj(1-zn)*lamba.real/ps

    print("Here is a small numerical check:\n")
    print("The numerical value of c_%i for z=%.2f%+.2fi:"%(n, zn.real, zn.imag))
    print(zeta.get_c(0))
    print("The value predicted by the formula:")
    print(lambda_x(zn.real, zn.imag) + 1j*lambda_y(zn.real, zn.imag))
    print("The value predicted by the simplified formula:")
    print(c_0_simplified)
