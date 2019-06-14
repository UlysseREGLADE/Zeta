#This code where used to studdy d_(x, y)(n) as decribed in the article

import numpy as np
import matplotlib.pyplot as plt
import sympy as sy
import poly_toolkite as po

def D(n, x, y):

    inv_sin = 1/np.sin(y*np.log((n+1)/n))
    inv_tan = 1/np.tan(y*np.log((n+2)/(n+1)))
    inv_n = 1/(1+1/n)**x

    return n*(y*inv_sin - y*inv_n*inv_tan-inv_n*(x-1))

def minoration_dl(n, x, y):

    inv_n = (4*x**3+3*x**2+4*x*y**2+7*y**2+1)/(12*n)

    return -inv_n + (x**2+y**2)/2

def limit(n, x, y):

    return np.ones((len(n)))*(x**2+y**2)/2

def D_eq(n, x, y):

        inv_sin = 1/po.sin(y*po.log_1(1/n, 3), 3)
        inv_tan = po.cos(y*po.log_1(1/(n+1), 3), 2)/po.sin(y*po.log_1(1/(n+1), 3), 3)
        inv_n = po.inv_1_(1/n, x, 3)

        return n*(y*inv_sin - y*inv_n*inv_tan-inv_n*(x-1))

def D_min(n, x, y):

        inv_sin = 1/po.sin(y*po.log_1(1/n, 3), 5)
        inv_tan = po.cos(y*po.log_1(1/(n+1), 3), 4)/po.sin(y*po.log_1(1/(n+1), 4), 3)
        inv_n = po.inv_1_(1/n, x, 4)

        return n*(y*inv_sin - y*inv_n*inv_tan-inv_n*(x-1))

def D_maj(n, x, y):

        inv_sin = 1/po.sin(y*po.log_1(1/n, 4), 3)
        inv_tan = po.cos(y*po.log_1(1/(n+1), 4), 2)/po.sin(y*po.log_1(1/(n+1), 3), 5)
        inv_n = po.inv_1_(1/n, x, 3)

        return n*(y*inv_sin - y*inv_n*inv_tan-inv_n*(x-1))

x, y = sy.symbols("x y", positive=True)
n = sy.symbols("n")

print("Working on D_min:")

exp_D_min = D_min(n, x, y) - (x**2+y**2)/2
print("assembeling fraction")
exp_D_min = sy.together(exp_D_min)
top, bot = sy.fraction(exp_D_min)
print("top degree: " + str(sy.degree(top, gen=n)))
print("bot degree: " + str(sy.degree(bot, gen=n)))

print('\n')

print("bot coefficients:")
poly_bot = sy.poly(bot, n)
for i in range(len(poly_bot.coeffs())):
    deg = sy.degree(bot, gen=n) - i
    coef = sy.factor(poly_bot.coeffs()[i])
    print("n**%i: "%(deg) + str(coef))

print('\n')

print("top coefficients:")
poly_top = sy.poly(top, n)
for i in range(len(poly_top.coeffs())):
    deg = sy.degree(top, gen=n) - i
    coef = sy.factor(poly_top.coeffs()[i])
    print("n**%i: "%(deg) + str(coef))

print('\n\n')

print("Working on D_maj:")

exp_D_maj = D_maj(n, x, y) - (x**2+y**2)/2
print("assembeling fraction")
exp_D_maj = sy.together(exp_D_maj)
top, bot = sy.fraction(exp_D_maj)
print("top degree: " + str(sy.degree(top, gen=n)))
print("bot degree: " + str(sy.degree(bot, gen=n)))

print('\n')

print("bot coefficients:")
poly_bot = sy.poly(bot, n)
for i in range(len(poly_bot.coeffs())):
    deg = sy.degree(bot, gen=n) - i
    coef = sy.factor(poly_bot.coeffs()[i])
    print("n**%i: "%(deg) + str(coef))

print('\n')

print("top coefficients:")
poly_top = sy.poly(top, n)
for i in range(len(poly_top.coeffs())):
    deg = sy.degree(top, gen=n) - i
    coef = sy.factor(poly_top.coeffs()[i])
    print("n**%i: "%(deg) + str(coef))


x=0.5
y=10

n = np.arange(10000)+1
plt.figure()
plt.title("plot of d_(%.2f, %.2f)"%(x, y))
plt.plot(D(n, x, y), label="d")
plt.plot(D_min(n, x, y), label="d_min")
plt.plot(D_maj(n, x, y), label="d_maj")
plt.plot(limit(n, x, y), label="lim")
plt.xlabel("n")
plt.ylabel("d_(%.2f, %.2f)(n)"%(x, y))
plt.legend()
plt.show(True)
