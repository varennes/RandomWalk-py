# Plot and create images of BPRW time evolution

import os
import numpy as np
import matplotlib as mp
import matplotlib.pyplot as plt

# input file with trajectory data
fIn = 'trj_6.dat'

with open( fIn, 'r') as f:
    x = []
    y = []
    i = 0
    for line in f:
        if (i % 2) == 0:
            x.append([float(v) for v in line.split()])
        else:
            y.append([float(v) for v in line.split()])
        i += 1

n = i / 2

tList = [ len(t) for t in x]
tmax  = max(tList)
tmin  = min(tList)

rmax = [ 0.0, 0.0]
rmin = [ min(x[0]), min(y[0])]
for i in range(len(x)):
    if max(x[i]) > rmax[0]:
        rmax[0] = max(x[i])
    if max(y[i]) > rmax[1]:
        rmax[1] = max(y[i])
    if min(x[i]) < rmin[0]:
        rmin[0] = min(x[i])
    if min(y[i]) < rmin[1]:
        rmin[1] = min(y[i])

# fig = plt.figure()
# for i in range(n):
#     plt.plot( x[i], y[i])
# plt.show()
# plt.close(fig)

for t in range(1,tmin):
    fig = plt.figure()
    for i in range(n):
        plt.plot( x[i][:(t+1)], y[i][:(t+1)], linewidth=2)

    plt.xlim( [rmin[0], rmax[0]])
    plt.ylim( [rmin[1], rmax[1]])
    plt.xticks([])
    plt.yticks([])
    plt.title('Biased Persistent Random Walk')

    fOut = 'trj_t_'
    if t < 10:
        fOut += '0'+str(t)
    else:
        fOut += str(t)
    fOut += '.png'
    plt.savefig(fOut, bbox='tight')
    plt.close(fig)
