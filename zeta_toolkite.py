#A little toolkite to easily manipulate the Riemann series

import numpy as np
from mpmath import zeta as true_zeta
import matplotlib.pyplot as plt

def intersec(z1, v1, z2, v2, n1=None):

    x1, y1 = z1.real, z1.imag
    x2, y2 = z2.real, z2.imag
    vx1, vy1 = v1.real, v1.imag
    vx2, vy2 = v2.real, v2.imag

    A = np.array([[vy1, -vx1],
                  [vy2, -vx2]])
    A_1 = np.array([[-vx2, vx1],
                    [-vy2, vy1]])
    B = np.array([[vy1*x1-vx1*y1],
                  [vy2*x2-vx2*y2]])

    print("Det A:")
    print(str(np.linalg.inv(A)))

    print(A_1*((((n1+1)*(n1+2))**1)/np.sin(-2*np.log((n1+2)/(n1+1)))))

    X = np.dot(np.linalg.inv(A), B)

    return X[0, 0] + 1j*X[1, 0]

class Zeta(object):

    def __init__(self, z):

        self.z = z

        self.zeta_n = [0, 1]

        self.zeta = true_zeta(self.z)

        if(self.z.imag>0):
            self.alpha = 3*np.pi/2 - np.arctan((1-self.z.real)/self.z.imag)
        elif(self.z.imag<0):
            self.alpha = np.pi/2 - np.arctan((1-self.z.real)/self.z.imag)

        self.gamma_n = []

        self.sum_arg_n = []

        self.theta_n = [0]

        self.c_n = []

    def get_zeta(self, n1=None):

        if(n1 is None):
            value = self.zeta
            return float(value.real)+float(value.imag)*1j

        n = n1

        if(type(n) is int):

            if(n < len(self.zeta_n)):
                return self.zeta_n[n]

            for i in range(len(self.zeta_n), n+2):
                self.zeta_n.append(self.zeta_n[-1] + 1/(i**self.z))

        else:

            if(np.max(n) < len(self.zeta_n)):
                return np.array(self.zeta_n)[n]

            for i in range(len(self.zeta_n), np.max(n)+2):
                self.zeta_n.append(self.zeta_n[-1] + 1/(i**self.z))

        return np.array(self.zeta_n)[n]

    def get_r(self, n):

        return 1/((n)**self.z)

    def approximate(self, n1):

        v1 = self.get_r(n1+1)*np.exp(1j*self.alpha)
        v2 = self.get_r(n1+2)*np.exp(1j*self.alpha)

        return intersec(self.get_zeta(n1), v1, self.get_zeta(n1+1), v2, n1=n1)

    def get_c(self, n1):


        n = n1

        if(type(n) is int):

            if(n < len(self.c_n)):
                return self.c_n[n]

            for i in range(len(self.c_n), n+2):
                self.c_n.append(self.approximate(i))

        else:

            if(np.max(n) < len(self.c_n)):
                return np.array(self.c_n)[n]

            for i in range(len(self.c_n), np.max(n)+2):
                self.c_n.append(self.approximate(i))

        return np.array(self.c_n)[n]

    def get_theta(self, n):

        if(type(n) is int):

            if(n < len(self.theta_n)):
                return self.theta_n[n]

            for i in range(len(self.theta_n), n+2):
                self.theta_n.append(np.angle(self.get_zeta(i+1)/self.get_zeta(i)))

        else:

            if(np.max(n) < len(self.theta_n)):
                return np.array(self.theta_n)[n]

            for i in range(len(self.theta_n), np.max(n)+2):
                self.theta_n.append(np.angle(self.get_zeta(i+1)/self.get_zeta(i)))

        return np.array(self.theta_n)[n]

        return np.angle(self.get_zeta(n+1)/self.get_zeta(n))

        return self.get_sum_arg(n+1) - self.get_sum_arg(n)

    def get_sum_arg(self, n):

        if(len(self.theta_n)-1 < n):

            self.get_theta(n)

        return np.sum(self.theta_n[:n+1])


    def get_gamma(self, n):

        return - self.z.imag*np.log(n+1) - self.get_sum_arg(n)

if(__name__=="__main__"):

    zeta = Zeta(1+2j)

    #print(zeta.get_zeta(10))
    #print(zeta.get_zeta(np.arange(10)+1))
    #print(len(zeta.zeta_n))
    print(zeta.get_c(1))
