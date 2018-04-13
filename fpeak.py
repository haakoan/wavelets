import numpy as np
import math
from matplotlib import pylab as plt
import itertools
from mpl_toolkits.axes_grid1 import make_axes_locatable


from matplotlib.ticker import AutoMinorLocator
from scipy.signal import filter_design as fd
plt.rc('text', usetex='True')
plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
plt.rcParams.update({'font.size': 36})#,'family':'monospace'})                                                                                                                                              
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['xtick.major.size'] = 10
plt.rcParams['xtick.major.width'] = 4
plt.rcParams['xtick.minor.size'] = 8
plt.rcParams['xtick.minor.width'] = 3

plt.rcParams['ytick.major.size'] = 10
plt.rcParams['ytick.major.width'] = 4
plt.rcParams['ytick.minor.size'] = 8
plt.rcParams['ytick.minor.width'] = 3

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True

plt.style.use('seaborn-colorblind')

#import random as rm
def fpeak(M,R,E):
    G = 6.67259e-8
    mn = 1.6749286e-24
    c = 2.99792458e10
    return 1/2/np.pi * G*M/R**2 * np.sqrt(1.1*mn/E) * (1 - G*M/R/c**2)**2
def fpeak2(args):
    return fpeak(*args)



msun = 1.99e33
kcm = 1e5
mev = 1/(624.15*1000)
Min = 1*msun
Rin = 100*kcm
Ein = 10*mev

xx = np.linspace(1*msun,2.5*msun,200)
yy = np.linspace(30*kcm,60*kcm,200)

# sj = []
# for x in xx:
#     kj = []
#     for y in yy:
#         kj.append(fpeak(x,y,Ein))
#     sj.append(kj)

#crd = list(itertools.product(*[xx,yy,[Ein]]))
#print(*crd)
#sj = list(map(fpeak2,crd))
f, ax = plt.subplots(nrows=2, ncols=2,figsize=(26,18))

v = np.linspace(50,1500,100)
for row in ax:
    for col in row:
        sj = []
        for x in xx:
            kj = []
            for y in yy:
                kj.append(fpeak(x,y,Ein))
            sj.append(kj)
        cl = col.contourf(xx/msun,yy/kcm,np.array(sj).T,v,cmap=plt.cm.plasma)
        CF = col.contour(xx/msun,yy/kcm,np.array(sj).T,10)
        Ein = Ein + 3*mev
        divider = make_axes_locatable(col)
        cax = divider.append_axes('right', size='7%', pad=0.13)
        col.set_xlabel(r'M [M$_\odot$]')
        col.set_ylabel('R [km]')
        tl = r'$<$E$_{\bar{\nu}_e}>$ = ' + str(math.ceil(Ein/mev)) + ' MeV'
        col.set_title(tl)
        cbar = f.colorbar(cl, cax=cax, orientation='vertical')
        cbar.add_lines(CF)
        col.tick_params(axis='x', pad=10)
        col.tick_params(axis='y', pad=10)
        col.yaxis.set_minor_locator(AutoMinorLocator(4))
        col.xaxis.set_minor_locator(AutoMinorLocator(5))
plt.subplots_adjust(hspace=0.4)
plt.subplots_adjust(wspace=0.4)

plt.savefig("spread.png")
print(fpeak(max(xx),min(yy),Ein))





# msun = 1.99e33
# kcm = 1e5
# mev = 1/(624.15*1000)
# Min = 1.2*msun
# Rin = 100*kcm
# Ein = 10*mev
# dm = np.random.uniform(0.3,0.3,5)*msun
# dr = np.random.uniform(-10,10,5)*kcm
# de = np.random.uniform(-4,4,5)*mev
# Ev = de + Ein
# Mv = dm + Min
# Rv = dr + Rin

# t = np.linspace(10,500,400)
# dt = t[1]-t[0]
# inital_conditions = list(itertools.product(*[Mv,Rv,Ev]))
# sj = list(map(fpeak2,inital_conditions))

# drdt = (30e-2)*dt*kcm
# dmdt = (0.075e-2)*dt*msun
# dedt = 0.012*dt*mev
# n = 0
# f = plt.figure()
# #plt.ylim(0,2000)
# while(n < 2):
#     gr = np.cumsum(-abs(np.random.normal(drdt, 15*drdt, len(t))))
#     gm = np.cumsum(abs(np.random.normal(dmdt, 15*dmdt, len(t))))
#     ge = np.cumsum(abs(np.random.normal(dedt, 15*dedt, len(t))))
#     print(gm[len(t)-1],gr[len(t)-1],ge[len(t)-1])
#     for x in inital_conditions:
#         vls = np.array((gm + x[0],gr + x[1],ge + x[2])).T
#         #print(vls.shape)
#         fs = list(map(fpeak2,vls))
#         #plt.plot(t,gm+x[0],'r.')
#         plt.plot(t,gr+x[1],'b.')
#         #plt.plot(t,ge+x[2],'k.')

#     n = n+1

# plt.savefig("spread.png")
