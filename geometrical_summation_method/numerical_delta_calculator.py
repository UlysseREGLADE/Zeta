#Code used to check numericaly the expression of delta_n(z)

import numpy as np
from zeta_toolkite import *
import matplotlib.pyplot as plt


x, y = 0.5, 2
z = x + y*1j
n=1
N=10

zeta = Zeta(z)

def delta_1(n, z):

    x, y = z.real, z.imag

    inv_sqrt = 1/np.sqrt(1+((1-z.real)/z.imag)**2)
    inv_sin = 1/np.sin(-z.imag*np.log((n+1)/n))
    inv_tan = 1/np.tan(-z.imag*np.log((n+2)/(n+1)))

    delta_s = inv_sqrt*(inv_sin*n**(-z.real) - (inv_tan+(1-z.real)/z.imag)*(n+1)**(-z.real))

    arg0 = zeta.alpha - y*np.log(n+1)

    return delta_s*np.exp(1j*arg0)

def delta_2(n, z):

    x, y = z.real, z.imag

    inv_sin = 1/np.sin(y*np.log((n+1)/n))
    inv_sin = 2*1j/(((n+1)/n)**(1j*y)-((n+1)/n)**(-1j*y))
    inv_tan = 1/np.tan(y*np.log((n+2)/(n+1)))
    inv_tan = 1j*(((n+2)/(n+1))**(1j*y)+((n+2)/(n+1))**(-1j*y))/(((n+2)/(n+1))**(1j*y)-((n+2)/(n+1))**(-1j*y))

    delta_s = inv_sin*n**(-x) - inv_tan*(n+1)**(-x) + ((1-x)/y)*(n+1)**(-x)
    arg0 = np.arctan((x-1)/y)
    exp = (1-np.conj(z))/np.abs(1-z)

    return (y/(np.abs(1-z)*(n+1)**(1j*y)))*delta_s*exp

def delta_3(n, z):

    x, y = z.real, z.imag

    inv_sin = 2*1j/(((n+1)/n)**(2j*y)-1)
    inv_tan = 1j*(((n+2)/(n+1))**(2j*y)+1)/(((n+2)/(n+1))**(2j*y)-1)

    delta_s = inv_sin*n**(-z) - inv_tan*(n+1)**(-z) + ((1-x)/y)*(n+1)**(-z)

    return (1-np.conj(z))*(y/(np.abs(1-z)**2))*delta_s

def sum_delta_1(N, z):

    x, y = z.real, z.imag

    zeta = Zeta(z)

    sum = zeta.get_zeta(N)-1
    mod = (1-np.conj(z))/(np.abs(1-z)**2)
    inv_sin = (y*2**(-1j*y))/np.sin(y*np.log(2))
    inv_tan = y/np.tan(y*np.log((N+2)/(N+1)))

    return sum + mod*(inv_sin - inv_tan*(N+1)**(-z) + (1-x)*(N+1)**(-z))

def c_1(N, z):

    x, y = z.real, z.imag

    zeta = Zeta(z)

    sum = zeta.get_zeta(N)
    mod = (1-np.conj(z))/(np.abs(1-z)**2)
    #inv_sin = (y*2**(-1j*y))/np.sin(y*np.log(2))
    inv_tan = y/np.tan(y*np.log((N+2)/(N+1)))

    return sum + mod*( - inv_tan*(N+1)**(-z) + (1-x)*(N+1)**(-z))

print("Numerical verification of value of delta_%i for z=%.2f%+.2fi:"%(n, x, y))
print("value from formula 1: " + str(delta_1(n, z)))
print("value from formula 2: " + str(delta_2(n, z)))
print("value from formula 3: " + str(delta_3(n, z)))
print("expected value: "+str(zeta.get_c(n)-zeta.get_c(n-1)))
print("\n")

print("Numercial check for the value of sum delta_1 to delta_%i  z=%.2f%+.2fi:"%(N, x, y))
print("value from formula 1: "+str(sum_delta_1(N, z)))
print("expected value: "+str(zeta.get_c(N)-zeta.get_c(0)))
print("\n")

print("Numercial check for the value of sum c_%i  z=%.2f%+.2fi:"%(N, x, y))
print("value from formula 1: "+str(c_1(N, z)))
print("expected value: "+str(zeta.get_c(N)))
print("\n")
