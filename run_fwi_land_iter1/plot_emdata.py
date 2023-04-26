import pandas as pd
import matplotlib.pyplot as plt
#plt.rcParams["font.family"] = "Helvetica"

import numpy as np
import cmath
import os

#-----------------------------------------------
iTx = []
iRx = []
ichrec = []
ifreq = []
dat = []

f =open('emf_0010.txt', 'r')
header = f.readline()
for line in f:
    line = line.strip()
    columns = line.split()
    iTx.append(int(columns[0]))
    iRx.append(int(columns[1]))
    ichrec.append(columns[2])
    ifreq.append(int(columns[3]))
    dat.append(float(columns[4]) +float(columns[5])*1j)
iTx = np.array(iTx)
iRx = np.array(iRx)
ichrec = np.array(ichrec)
ifreq = np.array(ifreq)
dat = np.array(dat)

amp = np.abs(dat) #nsrc*nfreq*nrec
pha = np.angle(dat, deg=True) #phase
 
iTx = []
iRx = []
ichrec = []
ifreq = []
dat = []

f =open('syn_0010.txt', 'r')
header = f.readline()
for line in f:
    line = line.strip()
    columns = line.split()
    iTx.append(int(columns[0]))
    iRx.append(int(columns[1]))
    ichrec.append(columns[2])
    ifreq.append(int(columns[3]))
    dat.append(float(columns[4]) +float(columns[5])*1j)
iTx = np.array(iTx)
iRx = np.array(iRx)
ichrec = np.array(ichrec)
ifreq = np.array(ifreq)
dat = np.array(dat)

amp2 = np.abs(dat) #nsrc*nfreq*nrec
pha2 = np.angle(dat, deg=True) #phase
 

plt.figure()
plt.subplot(211)
idx = (ifreq==1) & (ichrec=='Ex')
plt.plot(iRx[idx], amp[idx], 'r', label='obs Ex-0.25 Hz')
plt.plot(iRx[idx], amp2[idx], 'r--', label='syn Ex-0.25 Hz')
idx = (ifreq==2) & (ichrec=='Ex')
plt.plot(iRx[idx], amp[idx], 'g', label='obs Ex-0.75 Hz')
plt.plot(iRx[idx], amp2[idx], 'g--', label='syn Ex-0.75 Hz')

plt.title('(a)')
plt.xlabel('#iRx')
plt.yscale('log')
plt.ylabel('Amplitude (V/Am$^2$)')
plt.grid(color='k', linestyle='-')
plt.legend()


plt.subplot(212)
idx = (ifreq==1) & (ichrec=='Ex')
plt.plot(iRx[idx], pha[idx], 'r', label='obs Ex-0.25 Hz')
plt.plot(iRx[idx], pha2[idx], 'r--', label='syn Ex-0.25 Hz')
idx = (ifreq==2) & (ichrec=='Ex')
plt.plot(iRx[idx], pha[idx], 'g', label='obs Ex-0.75 Hz')
plt.plot(iRx[idx], pha2[idx], 'g--', label='syn Ex-0.75 Hz')


plt.title('(b)')
plt.xlabel('#iRx')
plt.ylabel('Degree ($^o$)')
plt.grid(color='k', linestyle='-')
plt.legend()

plt.tight_layout()
plt.savefig('comparison_obs_syn.pdf')
plt.show()

