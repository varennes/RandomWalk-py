# create and output trajectories of a BPRW

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

def getTheta( told, b, mu, k1, k2, h1, h2):
    m = 0
    while m < 1:
        y = 2.0*np.pi*( rd.random() - 0.5)
        u = rd.random()
        if u <= ( Ttheta2( y, told, b, mu, k1, k2, h1, h2)):
            m += 1
    return y


# BPRW parameters
b  = 0.1
mu = 0.0
k1 = 0.0
k2 = 5
h1 = 0.5
h2 = 0.5
thetaold = 0.0

l = 1.0
s = 3.75
tmax = 36.0
run = 3
dtSample = 2*l

fOut = 'trj_'+str(run)+'.dat'
print 'Trajectories output to %s' % fOut
with open( fOut, 'w') as f:
    pass

vrun = []
for i in xrange(run):
    x = [ 0.0]
    y = [ 0.0]
    r = []
    v = []
    t = 0.0
    tSample = 0.0 + dtSample
    r.append([ 0.0, 0.0, t])
    thetai = 2.0*np.pi*(rd.random()-0.5)
    # thetai = 0.0
    v.append([np.cos(thetai), np.sin(thetai), t])
    dt = np.random.exponential(l)
    t  = t + dt
    while t < tmax:
        r.append([ r[-1][j]+dt*v[-1][j]*s for j in range(2)])
        r[-1].append(t*s)
        thetaf = getTheta( thetai, b, mu, k1, k2, h1, h2)
        thetai = thetaf
        v.append([np.cos(thetaf), np.sin(thetaf), t])
        dt = np.random.exponential(l)
        t = t + dt
        if t >= tSample:
            x.append( r[-1][0])
            y.append( r[-1][1])
            tSample += dtSample
        # print 'run: '+str(i)+', t = '+str(t)
    dt = tmax - t
    r.append([ r[-1][j]+dt*v[-1][j] for j in range(2)])
    r[-1].append(tmax)

    vrun.append( v[-1][0:2])

    strX = ''
    strY = ''
    for i in range(len(x)):
        strX += str(x[i]) + ' '
        strY += str(y[i]) + ' '
    strX += '\n'
    strY += '\n'
    with open( fOut, 'a') as f:
        f.write(strX)
        f.write(strY)
