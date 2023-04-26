import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft,ifft


filename ="src_td"
f = open(filename, mode="rb")
dat1 = np.fromfile(f, dtype=np.float32)

filename ="src_new"
f = open(filename, mode="rb")
dat2 = np.fromfile(f, dtype=np.float32)

N = len(dat1)
x = [i*0.002 for i in range(N)]

yy=fft(dat1)
yf=abs(yy)  
xf = [i/(0.002*N) for i in range(N)]

plt.figure(figsize=(8,8))
plt.subplot(311)
plt.plot(x, dat1, 'k')
plt.xlabel('Time (s)')
plt.ylabel('S(t)')
plt.title('(a)')

plt.subplot(312)
plt.plot(xf[:10],yf[:10], 'k')
plt.xlabel('Frequency (Hz)')
plt.ylabel('S(w)')
plt.title('(b)')

plt.subplot(313)
plt.plot(x, dat2, 'k')
plt.xlabel('Time (s)')
plt.ylabel('Ss(t)')
plt.title('(c)')

plt.tight_layout()
plt.savefig('demo_1d.png')
plt.show()
