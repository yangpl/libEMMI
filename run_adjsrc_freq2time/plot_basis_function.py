import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['axes.xmargin'] = 0


f = open('basis_function', mode='rb')
dat = np.fromfile(f, dtype=np.float32)

t= [ i*0.004 for i in range(1000)]
traces = np.reshape(dat, (6,1000), order='C')

plt.figure(figsize=(12,9))

plt.subplot(311)
plt.plot(t, traces[0,:], 'k-', label='Real - b1')
plt.plot(t, traces[3,:], 'k--', label='Imag - b4')

plt.xlabel('Time (s)')
plt.legend()
plt.title('(a) 0.25 Hz')


plt.subplot(312)
plt.plot(t, traces[1,:], 'k-', label='Real - b2')
plt.plot(t, traces[4,:], 'k--', label='Imag - b5')

plt.xlabel('Time (s)')
plt.legend()
plt.title('(b) 0.75 Hz')


plt.subplot(313)
plt.plot(t, traces[2,:], 'k-', label='Real - b3')
plt.plot(t, traces[5,:], 'k--', label='Imag - b6')

plt.xlabel('Time (s)')
plt.legend()
plt.title('(c) 1.25 Hz')

plt.tight_layout(pad=1)
plt.savefig('basisfunction.png')
plt.show()
