import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
import pywt
from scipy.optimize import curve_fit
from scipy import signal
def wlt(t,A,f0,Q,t0,p):
    tau = Q/np.pi*0.5/f0
    return A*np.exp(-(t-t0)**2/tau**2)*np.cos(2.0*np.pi*f0*(t-t0) + p)


bh, ah = signal.butter(2, 0.2,'high')

data = np.loadtxt("processed_signal_s15a2o09_ls-plus.txt")
data = data[3478:3600]
#print(data)
x = np.linspace(0,len(data),len(data))
wavelet = pywt.ContinuousWavelet('morl')
psi, xx = wavelet.wavefun(level=10)
psi2 = signal.filtfilt(bh, ah, psi)
print(1/(xx[1]-xx[0]))
popt, pcov = curve_fit(wlt, xx, psi)



ftdata = wlt(xx, popt[0], popt[1], popt[2], popt[3], popt[4])
#print(popt)
#print((np.diag(pcov)))
#ft = (np.fft.fft(psi))
#ft2 = np.fft.fft(signal.filtfilt(bh, ah, ft))

f = plt.figure()
#plt.plot(x,data)
plt.plot(xx,psi)
plt.plot(xx,ftdata,'r--')
plt.savefig("plot.pdf")




# plt.plot(x,data)
# #plt.plot(xx-.8,coffs[0][0]*scaling*5e6)
# #plt.plot(xx-1.3,coffs[0][0]*scaling*2e9+coffs[1][0]*wavelet*1e6)
# #plt.plot(xx-2.7,-coffs[0][0]*scaling*4.2e12)
# plt.plot(xx-3.2,-scaling*9.2e7)

# #plt.plot(x, result.best_fit, 'r-')
# #plt.xlim(3400,3900)
# plt.savefig("plot.pdf")





# x = 10/len(x)*x
# w = pywt.Wavelet('morl')
# md = "zero"
# #cA, cD = pywt.dwt(data, wavelet=w,mode=md)
# #wl = pywt.idwt(cA, cD, wavelet=w, mode=md)
# #cA, cD = pywt.cwt(data, wavelet=w,mode=md)

# #coffs = pywt.wavedec(data, w, mode=md)
# #wl2 = pywt.waverec(coffs, w)
# coef, freqs=pywt.cwt(data,np.arange(1,129),'morl')
# scaling, wavelet, xx = w.wavefun()
# print(len(coffs[0]))
# print(coffs[0])
