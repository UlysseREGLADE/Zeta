#Code used to compute the assymptotic expansion of the Jacobian of delta_n(z)

import sympy as sy
import numpy as np
import pickle
import os

x = sy.symbols("x", real=True)
y = sy.symbols("y", real=True)
z = x + sy.I*y
n = sy.symbols("n")

def delta_x(n, z):

    x, y = sy.re(z), sy.im(z)

    inv_sqrt = (1/sy.sqrt(1+((1-x)/y)**2))
    inv_sin = (n**-x/sy.sin(-y*sy.log((n+1)/n)))
    inv_tan = ((n+1)**-x/sy.tan(-y*sy.log((n+2)/(n+1))))
    const = ((1-x)*(n+1)**-x/y)
    cos = sy.cos(sy.pi/2-sy.atan((1-x)/y)-y*sy.log(n+1))

    #cos = inv_sqrt*(sy.cos(y*sy.log(n+1))*(1-x)/y + sy.sin(y*sy.log(n+1)))

    return inv_sqrt*(inv_sin-inv_tan-const)*cos

def delta_y(n, z):

    x, y = sy.re(z), sy.im(z)

    inv_sqrt = (1/sy.sqrt(1+((1-x)/y)**2))
    inv_sin = (n**-x/sy.sin(-y*sy.log((n+1)/n)))
    inv_tan = ((n+1)**-x/sy.tan(-y*sy.log((n+2)/(n+1))))
    const = ((1-x)*(n+1)**-x/y)
    sin = sy.sin(sy.pi/2-sy.atan((1-x)/y)-y*sy.log(n+1))

    #sin = inv_sqrt*(sy.cos(y*sy.log(n+1)) - sy.sin(y*sy.log(n+1))*(1-x)/y)

    return inv_sqrt*(inv_sin-inv_tan-const)*sin

def delta_r(n, z):

    x, y = sy.re(z), sy.im(z)

    inv_sqrt = (1/sy.sqrt(1+((1-x)/y)**2))
    inv_sin = (n**-x/sy.sin(-y*sy.log((n+1)/n)))
    inv_tan = ((n+1)**-x/sy.tan(-y*sy.log((n+2)/(n+1))))
    const = ((1-x)*(n+1)**-x/y)

    return sy.Abs(inv_sqrt*(inv_sin-inv_tan-const))

def delta_theta(n, z):

    x, y = sy.re(z), sy.im(z)

    return sy.pi/2 - sy.atan((1-x)/y) - y*sy.log(n+1)

delta_x = ((1/sy.sqrt(1+((1-x)/y)**2))*((n**-x/sy.sin(-y*sy.log((n+1)/n)))-((n+1)**-x/sy.tan(-y*sy.log((n+2)/(n+1))))-((1-x)*(n+1)**-x/y)))*sy.cos(sy.pi/2-sy.atan((1-x)/y)-y*sy.log(n+1))
delta_y = ((1/sy.sqrt(1+((1-x)/y)**2))*((n**-x/sy.sin(-y*sy.log((n+1)/n)))-((n+1)**-x/sy.tan(-y*sy.log((n+2)/(n+1))))-((1-x)*(n+1)**-x/y)))*sy.sin(sy.pi/2-sy.atan((1-x)/y)-y*sy.log(n+1))
delta_r = (1/sy.sqrt(1+((1-x)/y)**2))*((1/sy.sin(-y*sy.log((n+1)/n)))-((1+1/n)**-x/sy.tan(-y*sy.log((n+2)/(n+1))))-((1-x)*(1+1/n)**-x/y))

J_r_y = sy.diff(delta_r, y)
J_r_x = sy.diff(delta_r, x)

print("J_r_x")
sy.pprint(J_r_x)
print(sy.latex(J_r_x) + '\n')
print("J_r_y")
sy.pprint(J_r_y)
print(sy.latex(J_r_y) + '\n')

dlx = J_r_x.series(n, sy.oo, 2)
dly = J_r_y.series(n, sy.oo, 2)

print("dlx")
sy.pprint(dlx)
print("dly")
sy.pprint(dly)

dlx_O = sy.simplify(n*dlx.removeO())
print("dlx_O")
sy.pprint(dlx_O)
dly_O = sy.simplify(n*dly.removeO())
print("dly_O")
sy.pprint(dly_O)
