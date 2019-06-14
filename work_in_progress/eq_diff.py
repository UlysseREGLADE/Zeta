import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
from mpmath import zeta

M, N = 200, 200
mp.mp.prec = 1000

def r(n, m):

    return (-n)/(1-z-m)
x = mp.mpc(0.5)
y = mp.mpc('14.1347251417346937904572519835624702707842571156992431756855674601499634298092567649490103931715610127792')
z = x + 1j*y
print(z)

print(np.angle(np.array(zeta(z)).astype(complex)))

x, y = z.real, z.imag

zeta_m_n = []

for m in range(M):

    zeta_n  = [mp.mpc(0)]

    if(m==0):

        for n in range(1, N):

            zeta_n.append(zeta_n[-1] + mp.mpc(n**-z))

    else:

        for n in range(1, len(zeta_m_n[-1])-1):

            zeta=zeta_m_n[m-1][n]+(zeta_m_n[m-1][n+1]-zeta_m_n[m-1][n])*r(n,m-1)

            zeta_n.append(zeta)

    zeta_m_n.append(np.array(zeta_n))

zeta_m_n = np.array(zeta_m_n)

Z = np.zeros((M, N), dtype=complex)
X, Y = np.arange(N), np.arange(M)
X, Y = np.meshgrid(X, Y)

plt.figure()
plt.title("zeta_n^m(z) for z = "+str(z))
for m in range(len(zeta_m_n)):
    zeta_m_n[m] = zeta_m_n[m].astype(complex)
    Z[m, :len(zeta_m_n[m])] = zeta_m_n[m]
    if(m<6):
        plt.scatter(zeta_m_n[m].real,
                    zeta_m_n[m].imag,
                    s=0.3,
                    label="zeta_n^%i"%(m))
plt.plot(mp.zeta(z).real, mp.zeta(z).imag, '+', label="zeta(z)")
plt.axis('equal')
plt.legend()

def format_coord_angle(x, y):
    try:
        int_x, int_y = int(x), int(y)
        z = np.angle(Z[int_y, int_x])
        return "x: {}, y: {}, angle: {}".format(int_x,int_y,z)
    except:
        return "x: {}, y: {}".format(x,y)

def format_coord_abs(x, y):
    try:
        int_x, int_y = int(x), int(y)
        z = np.abs(Z[int_y, int_x])
        return "x: {}, y: {}, abs: {}".format(int_x,int_y,z)
    except:
        return "x: {}, y: {}".format(x,y)


f, (ax1, ax2) = plt.subplots(1, 2)
cmap = plt.cm.hsv
ax1.pcolormesh(X, Y, np.angle(Z), cmap=cmap)
ax1.format_coord = format_coord_angle
#ax1.set_aspect('equal')
ax2.pcolormesh(X, Y, np.log(np.abs(Z)))
ax2.format_coord = format_coord_abs
#ax2.set_aspect('equal')
plt.show()
