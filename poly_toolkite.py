#A little toolkite to facilitate the calculation of assymptotic expansions

import sympy as sy
import numpy as np
import matplotlib.pyplot as plt

def fact(i):

    ret = 1

    for i in range(1, i+1):

        ret *= i

    return ret


def sin(x, n):

    i = n//2 + (n%2)

    ret = 0

    for j in range(i):

        ret += (-1)**j*x**(2*j+1)/fact(2*j+1)

    return ret

def cos(x, n):

    i = n//2

    ret = 1

    for j in range(i):

        ret += (-1)**(j+1)*x**(2*(j+1))/fact(2*(j+1))

    return ret

def log_1(x, n):

    ret = 0

    for i in range(n):

        ret += (-1)**i*x**(i+1)/(i+1)

    return ret

def inv_1_(x, alpha, n):

    ret = 1

    ALPHA = -alpha

    for i in range(n):

        ret += x**(i+1)*ALPHA/fact(i+1)
        ALPHA *= (-alpha-i-1)

    return ret

if(__name__ == "__main__"):

    x = sy.symbols('x')
    n=4
    alpha = 2

    sy.pprint(inv_1_(x, alpha, n))

    x = np.linspace(0, 2, 100)

    plt.figure()
    plt.plot(inv_1_(x, alpha, n))
    plt.plot(1/(1+x)**alpha)
    plt.show()
