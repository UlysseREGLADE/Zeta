import sympy as sy
import mpmath as mp
import poly_toolkite as po

M_max = 6

M_max += 1

N, z = sy.symbols('N z')

R = [sy.S('0')]

for M in range(M_max):

    new_R = R[-1]
    new_R += (R[-1].subs(N, N+1)-R[-1]+(N+1)**(-z))*N/(M+z-1)

    R.append(new_R)

x = sy.symbols('x')

exp_approxs = []

for m in range(M_max):

    for i in range(M_max):

        I = i+1

        R[m] = R[m].subs((N+I)**-z, x**I)

    R[m] = sy.apart(R[m], x)

    print("Corrective terme of order: M=%i\n"%(m))

    coeffs = sy.poly(R[m], x).coeffs()

    exp_approx = 1

    for i in range(len(coeffs)):

        I = m-i

        coeff = coeffs[i]
        coeff = coeff.subs(N, 1)
        exp = po.exp(-z*sy.log(I+1), m-1)
        exp_approx += coeff*exp

        sy.pprint(1/(N+I)**z)
        #sy.pprint(po.exp(-z*sy.log(I+1), m-1))
        sy.pprint(coeff)
        print('\n')

    exp_approxs.append(exp_approx)

    print('\n')

frac = sy.together(exp_approxs[-1])
top, bot = sy.fraction(frac)

poly_top = sy.poly(top, z)
sy.pprint(poly_top)
sy.pprint(sy.factor(poly_top))
sy.pprint(poly_top.coeffs())

fact = []
for i in range(len(poly_top.coeffs())):

    fact.append(sy.N(poly_top.coeffs()[i]))
print(fact)
# roots = mp.polyroots(fact, error=True)
# for root in roots:
#     print(root)
