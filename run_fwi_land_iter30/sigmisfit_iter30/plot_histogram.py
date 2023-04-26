import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob

plt.rcParams['axes.xmargin'] = 0

iTx = np.zeros((1,)) #0 rows, 1 column empty np array
iRx = np.zeros((1,)) #0 rows, 1 column empty np array
ichrec = []
ifreq = np.zeros((1,)) #0 rows, 1 column empty np array
obj = np.zeros((1,)) #0 rows, 1 column empty np array

filelist = glob.glob('sig_misfit_00*.txt')
for fname in filelist:
    print(fname)
    f =open(fname, 'r')
    header = f.readline()
    for line in f:
        line = line.strip()
        columns = line.split()

        iTx = np.hstack((iTx, int(columns[0])))
        iRx = np.hstack((iRx, int(columns[1])))
        ichrec.append(columns[2])
        ifreq = np.hstack((ifreq, int(columns[3])))
        obj = np.hstack((obj, float(columns[4])))
ichrec = np.array(ichrec)


df = pd.DataFrame(data=obj)
df.plot(kind='hist',bins=100, range=[0,30], grid=True, legend=False, orientation='vertical', color='gray', figsize=(8,4))


plt.xlabel('$|W(d_{obs}-d_{syn})| $')
plt.title('(b)')

plt.tight_layout(pad=0)
plt.savefig('histogram.eps')
plt.show()
