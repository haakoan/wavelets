import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
import pywt
from scipy.optimize import curve_fit
data = np.loadtxt("processed_signal_s15a2o09_ls-plus.txt")
data = data[3478:3600]
#print(data)
x = np.linspace(0,len(data),len(data))
x = 10/len(x)*x
w = pywt.Wavelet('db15')
md = "zero"
cA, cD = pywt.dwt(data, wavelet=w,mode=md)
wl = pywt.idwt(cA, cD, wavelet=w, mode=md)

coffs = pywt.wavedec(data, w, mode=md, level=2)
wl2 = pywt.waverec(coffs, w)

scaling, wavelet, xx = w.wavefun()
print(len(coffs[0]))
print(coffs[0])


print(xx)
f = plt.figure()
plt.plot(x,data)
#plt.plot(xx-.8,coffs[0][0]*scaling*5e6)
#plt.plot(xx-1.3,coffs[0][0]*scaling*2e9+coffs[1][0]*wavelet*1e6)
plt.plot(xx-2.7,-coffs[0][0]*scaling*4.2e12)

#plt.plot(x, result.best_fit, 'r-')
#plt.xlim(3400,3900)
plt.savefig("plot.pdf")
