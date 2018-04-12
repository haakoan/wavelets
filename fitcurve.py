import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
import pywt
from scipy.optimize import curve_fit
def gaussian(x, amp, cen, wid):
    return amp * np.sin(cen*x)*np.exp(-x)

def har(x, a, d, f, p, k, a2, l, l2):
    return a*np.exp(d*x)*np.cos(f*x + p) + k + a2*np.exp((x-l)**2*l2)


#data = np.loadtxt("processed_signal_s15a2o05_ls-plus.txt")
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
# fit = np.polyfit(x,data,40)
# print(len(data))
# print(len(x))
# p = np.poly1d(fit)


# gmodel = Model(har)
# result = gmodel.fit(data, x=x, a=5, d=-1, f=20, p = 0, k = 1e-20,a2=0.1e-20,l=100,l2=-0.1)
# print(result.fit_report())

#def f_fit(t,w,cof,N):
#    funs
#    for l in range(0,N):
#        for i in range(0,2**l):
#        funs = funs

print(xx)
f = plt.figure()
plt.plot(x,data)
#plt.plot(xx-.8,coffs[0][0]*scaling*5e6)
#plt.plot(xx-1.3,coffs[0][0]*scaling*2e9+coffs[1][0]*wavelet*1e6)
plt.plot(xx-2.7,-coffs[0][0]*scaling*4.2e12)

#plt.plot(x, result.best_fit, 'r-')
#plt.xlim(3400,3900)
plt.savefig("plot.pdf")
