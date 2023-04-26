import pandas as pd
import matplotlib.pyplot as plt


plt.rcParams['axes.xmargin'] = 0


import numpy as np
import cmath
import os

#-----------------------------------------------
iteration = []
obj = []

f =open('conv_history.txt', 'r')
header = f.readline()
for line in f:
    line = line.strip()
    columns = line.split()
    #print(columns)
    iteration.append(int(columns[1]))
    obj.append(float(columns[4]))
iteration = np.array(iteration)
obj = np.array(obj)

#plt.grid(True, which="both", ls="-")
#plt.xticks(range(0,21,5))
plt.yscale('log')
plt.xlabel('# CG iteration')
plt.ylabel('$f_k/f_1$')
plt.plot(iteration, obj, 'k')
plt.savefig('CG_conv.png')
plt.show()
