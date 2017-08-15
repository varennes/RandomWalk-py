import matplotlib.pyplot as plt
from scipy import special
import numpy as np
import random as rd
np.random.seed()
rd.seed()

def hphi2( phi, mu, k1, k2, h1, h2):
    h1top = np.exp(k1*np.cos(phi-mu))
    h1bot = 2.0*np.pi*special.iv(0, k1)
    h2top = np.exp(k2*np.cos(phi-mu))
    h2bot = 2.0*np.pi*special.iv(0, k2)
    h     = h1*(h1top/h1bot) + h2*(h2top/h2bot)
    return h

def ktheta( theta, b):
    return b*np.cos(theta)

# reorientation kernel
def Ttheta2( tnew, told, b, mu, k1, k2, h1, h2):
    phi = tnew - told
    h   = hphi2( phi, mu, k1, k2, h1, h2)
    k   = ktheta( tnew, b)
    return h+k

b  = 0.1
mu = 0.0
k1 = 0.0
k2 = 5
h1 = 0.75
h2 = 0.25
thetaold = 0.0

n  = 100
xl = [-np.pi + 2.0*np.pi*float(i)/float(n) for i in range(n+1)]

yl = [ Ttheta2( x, thetaold, b, mu, k1, k2, h1, h2) for x in xl]
yh = [ hphi2( x-thetaold, mu, k1, k2, h1, h2) for x in xl]
yk = [ ktheta( x, b) for x in xl]

if min(yl) < 0.0:
    print 'ERROR: chosen b value too large, T < 0'
    b = round(min(yh),4)
    print '       new b value = '+str(b)

    yl = [ Ttheta2( x, thetaold, b, mu, k1, k2, h1, h2) for x in xl]
    yh = [ hphi2( x-thetaold, mu, k1, k2, h1, h2) for x in xl]
    yk = [ ktheta( x, b) for x in xl]


m = 0
M = 10000
sample = []
while m < M:
    y = 2.0*np.pi*( rd.random() - 0.5)
    u = rd.random()
    if u <= ( Ttheta2( y, thetaold, b, mu, k1, k2, h1, h2)):
        sample.append(y)
        m += 1

numbins = 50
nl, binl = np.histogram( sample, bins=numbins)
wl = binl[1]-binl[0]
norm = sum( nl*np.diff(binl) )
nl = nl / norm
nl = nl

plt.figure
# lT = r"$T_2(\theta,\theta')$"
# lh = r"$h_2(\theta-\theta')$"
lT = r"$T(\theta,\theta')$"
lh = r"$h(\theta-\theta')$"
lk = r'$k(\theta)$'
plt.bar( binl[:-1], nl, width=wl, color='lightgrey')
plt.plot( xl, yl, '-r',  linewidth=3, label=lT)
plt.plot( xl, yh, '-b', linewidth=2, label=lh)
plt.plot( xl, yk, '-g', linewidth=2, label=lk)
plt.legend(loc='best')

plt.xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi], [r'$-\pi$', r'$-\pi/2$', '0', r'$\pi/2$', r'$\pi$'], fontsize=11)

# ltitle =  r"$\theta'=0.0$, $\mu=0.0$, $\kappa_1=$ "+str(k1)+r", $\kappa_2=$ "+str(k2)
# ltitle += r", $h_1$ = "+str(h1)+r", $h_2$ = "+str(h2)+", b = "+str(b)

ltitle =  "Reorientation Probability Density"

plt.title(ltitle)
plt.savefig('ReorientProb_1.png', bbox='tight')
plt.show()
