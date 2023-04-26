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
obj = []
x = []
y = []
z = []


f =open('sig_misfit_0010.txt', 'r')
header = f.readline()
for line in f:
    line = line.strip()
    columns = line.split()
    iTx.append(int(columns[0]))
    iRx.append(int(columns[1]))
    ichrec.append(columns[2])
    ifreq.append(int(columns[3]))
    obj.append(float(columns[4]))
    x.append(float(columns[5]))
    y.append(float(columns[6]))
    z.append(float(columns[7]))
iTx = np.array(iTx)
iRx = np.array(iRx)
ichrec = np.array(ichrec)
ifreq = np.array(ifreq)
obj = np.array(obj)
x = np.array(x)
y = np.array(y)
z = np.array(z)

plt.figure(figsize=(12,9))
plt.subplot(221)
idx = (ifreq==1) & (ichrec=='Ex')
plt.scatter(x[idx], y[idx], s=50, c=obj[idx], cmap='rainbow')
plt.xlabel('X(m)')
plt.ylabel('Y(m)')
plt.title('ifreq=1, Ex')
plt.colorbar(location='right')


plt.subplot(222)
idx = (ifreq==2) & (ichrec=='Ex')
plt.scatter(x[idx], y[idx], s=50, c=obj[idx], cmap='rainbow')
plt.xlabel('X(m)')
plt.ylabel('Y(m)')
plt.title('ifreq=2, Ex')
plt.colorbar(location='right')



plt.subplot(223)
idx = (ifreq==1) & (ichrec=='Ey')
plt.scatter(x[idx], y[idx], s=50, c=obj[idx], cmap='rainbow')
plt.xlabel('X(m)')
plt.ylabel('Y(m)')
plt.title('ifreq=1, Ey')
plt.colorbar(location='right')


plt.subplot(224)
idx = (ifreq==2) & (ichrec=='Ey')
plt.scatter(x[idx], y[idx], s=50, c=obj[idx], cmap='rainbow')
plt.xlabel('X(m)')
plt.ylabel('Y(m)')
plt.title('ifreq=2, Ey')
plt.colorbar(location='right')

plt.tight_layout(pad=0)
plt.savefig('sigmisfit_scatter.png', bbox_inches='tight')
plt.show()


plt.figure()
idx = (ifreq==1) & (ichrec=='Ex')
yy = obj[idx]
idx = (ifreq==2) & (ichrec=='Ex')
yy += obj[idx]
idx = (ifreq==1) & (ichrec=='Ey')
yy += obj[idx]
idx = (ifreq==2) & (ichrec=='Ey')
yy += obj[idx]

plt.scatter(x[idx], y[idx], s=50, c=yy, cmap='rainbow')
plt.xlabel('X(m)')
plt.ylabel('Y(m)')
plt.title('$\sum obj$')
plt.colorbar()
plt.show()
