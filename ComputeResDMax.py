#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
#reload(sys)
#sys.setdefaultencoding('utf8')
# /---------------------------------|
'''
This script compute maximum euclidean distance of each residue and generate
a plot. In addition, a Haar wavelet transformation is provided to facilitate
the identification of most flexible residues.
'''
import math
from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons
import numpy as np

# Define Plot Language
lang = "EN"
#lang = "PT"
sign = lambda x: math.copysign(1, x)

# Defining Functions ##

# Peak detection test


def peak(z, i, t):
    '''
    Decide if value z(i) is a peak based on some conditions
    z = value
    i = indice of value z
    t = threshold
    '''

    if ((z[i] >= z[i-1]) or (z[i] >= z[i+1])) and (z[i] >= t):
        return True
    else:
        return False


# Haar wavelet transform
def haar_1d(z):
    '''    '''
    n = len(z)
    s = math.sqrt(2.0)
    x = np.array(z)
    y = np.zeros(n)

    k = 1
    while(k*2 <= n):
        k = k * 2

    while(1 < k):
        k = int(k/2)
        for i in range(k):
            y[i] = (x[2*i] + x[2*i+1])/s
            y[i+k] = (x[2*i] - x[2*i+1])/s

        for i in range(k*2):
            x[i] = y[i]

    return(x)


# Inverse Haar wavelet transform
def haar_1d_inverse(z):
    ''' '''
    n = len(z)
    s = math.sqrt(2.0)
    x = np.array(z)
    y = np.zeros(n)

    k = 1
    while(k*2 <= n):
        for i in range(k):
            y[2*i] = (x[i] + x[i+k])/s
            y[2*i+1] = (x[i] - x[i+k])/s

        for i in range(k*2):
            x[i] = y[i]

        k = k*2
    return(x)


########################################
# creating ssv file
in_file = open('DiffGeoSpec.ssv', 'r')
lines = in_file.readlines()
in_file.close()

# variables to store:
pos = []    # residues position
max_ct = []  # max curvature
max_tr = []  # max torsion
max_wt = []  # max withing
max_al = []  # max arc length
max_dt = []  # max distance

for line in lines:
    tokens = line.split()
    tokens = list(map(float, tokens))
    pos.append(tokens[0])
    max_ct.append(tokens[1])
    max_tr.append(tokens[2])
    max_wt.append(tokens[3])
    max_al.append(tokens[4])

    # euclidean distance

    # dt = sqrt(tokens[1]**2 + tokens[2]**2 + tokens[3]**2 + tokens[4]**2)
    dt = sqrt(tokens[1]**2 + tokens[2]**2)
    max_dt.append(dt)

data = max_dt

mu = np.mean(data)  # data mean
sigma = np.std(data)   # data standard deviation

# write a z-score csv
maxdt_outf = open("maxdt.csv", "w")

for n, res_n in enumerate(pos):
    maxdt_outf.write(str(res_n)+","+str(data[n])+"\n")
maxdt_outf.close()

# z
z = []
for i in range(len(pos)):
    test = (data[i]-mu)/sigma
    if data[i] == 0.0:
        test = 0.0
    z.append(test)

# write a z-score csv
z_outf = open("z-score.csv", "w")

for n, res_n in enumerate(pos):
    z_outf.write(str(res_n)+","+str(z[n])+"\n")
z_outf.close()

# Haar wavelet data analysis
limit = 0.4
thres = 0.4
zt = haar_1d(z)
zt = list(map(lambda v: v if (abs(v) > limit) else 0.0, zt)) # High-pass filter
zi = haar_1d_inverse(zt)

# Test zi for peaks
pos_p = [] #residue number
value = [] # values
peaks = []
for i in range(1, len(zt)-1):
    if (peak(zi, i, thres)):
            pos_p.append(pos[i])
            value.append(z[i])
            peaks.append(zi[i])

# writing residues list

print("@> writing aa_list.csv (limit =", limit, "; threshold =", thres, ")")
print("@ ", len(pos_p), "were found")
print("@ ", pos_p)
aa_list_file = open('aa_list.csv', 'w')
#aa_list_file.write('#list of representative position to consider')

for p in pos_p:
    aa_list_file.write(str(int(p))+"\n")


# Plotting
ax1 = subplot(211)
k, = plot(pos, zi)
j, = plot(pos_p, peaks, 'r.')
ylabel('Haar transform')
grid(True)
ax2 = subplot(212)
subplots_adjust(bottom=0.25)
plot(pos, z)
l, = plot(pos_p, value, 'r.')
grid(True)
if lang is "EN":
    xlabel('Residue Position')
    ylabel(r'$\kappa$/$\tau$ max distance (Z-Score)')
if lang is "PT":
    xlabel('Resíduo')
    ylabel(r'$\kappa$/$\tau$ distância máxima (Z-Score)')
#axis([0, 1, -10, 10])

axcolor = 'lightgoldenrodyellow'
axlimit = axes([0.25, 0.10, 0.65, 0.03])#, axisbg=axcolor)
axthres = axes([0.25, 0.15, 0.65, 0.03])#, axisbg=axcolor)

slimit = Slider(axlimit, 'Limit', 0.1, 4.0, valinit=0.4)
sthres = Slider(axthres, 'Thres', 0.1, 2.0, valinit=0.4)

def update(val):
    # Haar wavelet data analysis
    limit = slimit.val
    thres = sthres.val
    zt = haar_1d(z)
    zt = list(map(lambda v: v if (abs(v) > limit) else 0.0, zt)) # High-pass filter
    zi = haar_1d_inverse(zt)

    # Test zi for peaks
    pos_p = []
    value = []
    peaks = []
    for i in range(1, len(zt)-1):
        if (peak(zi, i, thres)):
            pos_p.append(pos[i])
            value.append(z[i])
            peaks.append(zi[i])
    k.set_ydata(zi)
    j.set_xdata(pos_p)
    j.set_ydata(peaks)
    l.set_xdata(pos_p)
    l.set_ydata(value)
    draw()

sthres.on_changed(update)
slimit.on_changed(update)

resetax = axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
def reset(event):
    sthres.reset()
    slimit.reset()
button.on_clicked(reset)
show()
