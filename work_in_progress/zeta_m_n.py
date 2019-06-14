import sympy as sy

M_max = 6

M_max += 1

N, z = sy.symbols('N z')

R = [sy.S('0')]

for M in range(M_max):

    new_R = R[-1]
    new_R += (R[-1].subs(N, N+1)-R[-1]+(N+1)**(-z))*N/(M+z-1)

    R.append(new_R)

x = sy.symbols('x')

for m in range(M_max):

    for i in range(M_max):

        I = i+1

        R[m] = R[m].subs((N+I)**-z, x**I)

    R[m] = sy.apart(R[m], x)

    print("Corrective terme of order: M=%i\n"%(m))

    coeffs = sy.poly(R[m], x).coeffs()

    for i in range(len(coeffs)):

        I = m-i

        coeff = coeffs[i]

        sy.pprint(1/(N+I)**z)
        if(coeff != 0):
            if(m!=1):
                aparts = sy.apart(coeff, z).args
                for apart in aparts:
                    sy.pprint(apart.subs(N, 1))
            else:
                apart = sy.apart(coeff, z)
                sy.pprint(apart.subs(N, 1))
        else:
            sy.pprint(0)
        print('\n')

    print('\n')
