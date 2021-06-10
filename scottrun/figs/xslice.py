import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os
from pylab import *
from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
sformatter=ScalarFormatter(useOffset=True,useMathText=True)
sformatter.set_scientific(True)
sformatter.set_powerlimits((-2,3))

Npanels=3

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 12}
plt.rc('font', **font)
plt.rc('text', usetex=False)
linestyles = [
  '-', '--', 'dotted', 'dash-dotted'
]
markerstyles = [
  u's',
  u'o',
  u'^',
  r'$\\circlearrowleft$',
  r'$\\clubsuit$',
  r'$\\checkmark$'
]
colors = ('b', 'g', 'r', 'c', 'm', 'y', 'k')

plt.figure(figsize=(6,10))
fig = plt.figure(1)

xx0=0.18
yy0=0.07
ww0=1.0-xx0-0.02
hh0=1.0-yy0-0.02

c2=5.0/3.0
Lmax=100.0
k=2.0*pi/Lmax
omega0=sqrt(c2)*k
print('omega=',omega0,' period=',2.0*pi/omega0)

timearray=['0','100','200','300','400','500','600','700','800','900','1000']

modeldir='../../output_posterior_corrected/default/'

ipanel=0
ymax0=1.0012
ymin0=0.9988
dely=0.0004
ax = fig.add_axes([xx0,yy0+ipanel*hh0/Npanels,ww0,hh0/Npanels])
ax.set_xticks(np.arange(0,1.2*Lmax,0.2*Lmax), minor=False)
ax.set_xticks(np.arange(0,1.2*Lmax,0.1*Lmax), minor=True)
ax.set_xticklabels(np.arange(0,1.2*Lmax,0.2*Lmax), minor=False, size=14)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
ax.set_yticks(np.arange(ymin0,ymax0,dely), minor=False)
ax.set_yticklabels(np.arange(ymin0,ymax0,dely), minor=False, family='serif', size=14)
ax.set_yticks(np.arange(ymin0,ymax0,0.5*dely), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.4f'))
plt.ylabel("$\\rho$")
plt.xlabel('$x$',fontsize=24)
plt.ylim(ymin0,ymax0)
plt.xlim(0.0,Lmax)

for time in timearray:
  tt=float(time)
  filename='../output/xslice_t'+str(time)+'.dat'
  print('filename=',filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  yguess=arange(0,len(x),1.0)
  rhoB=data[1]
  jB=data[2]
  epsilon=data[3]
  ux=data[4]
  pressure=data[5]
  plt.plot(x,rhoB,linestyle='-',linewidth=3,color='r')
  #plt.clf()
  for i in range(0,len(x)):
    yguess[i]=1.0+0.001*cos(omega0*tt)*cos(k*x[i])
  plt.plot(x,yguess,linestyle='--',linewidth=3,color='k')
  
ipanel=1
ymin1=0.996
ymax1=1.004
dely=0.002
ax = fig.add_axes([xx0,yy0+ipanel*hh0/Npanels,ww0,hh0/Npanels])
ax.set_xticks(np.arange(0,1.2*Lmax,0.2*Lmax), minor=False)
ax.set_xticks(np.arange(0,1.2*Lmax,0.1*Lmax), minor=True)
ax.set_xticklabels([], minor=False)
ax.set_yticks(np.arange(ymin1,ymax1,dely), minor=False)
ax.set_yticklabels(np.arange(ymin1,ymax1,dely), minor=False, family='serif', size=14)
ax.set_yticks(np.arange(ymin1,ymax1,0.5*dely), minor=True)
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
plt.ylabel('$P$')
plt.ylim(ymin1,ymax1)
plt.xlim(0.0,Lmax)

for time in timearray:
  tt=float(time)
  filename='../output/xslice_t'+str(time)+'.dat'
  print('filename=',filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  rhoB=data[1]
  jB=data[2]
  epsilon=data[3]
  ux=data[4]
  pressure=data[5]
  ax.plot(x,pressure,linestyle='-',linewidth=3,color='g')
  for i in range(0,len(x)):
    yguess[i]=(1.0+0.001*cos(omega0*tt)*cos(k*x[i]))**(1.666666667)
  plt.plot(x,yguess,linestyle='--',linewidth=3,color='k')

ipanel=2
ymax2=0.002
ymin2=-ymax2
dely=0.001
ax = fig.add_axes([xx0,yy0+ipanel*hh0/Npanels,ww0,hh0/Npanels])
ax.set_xticks(np.arange(0,1.2*Lmax,0.2*Lmax), minor=False)
ax.set_xticks(np.arange(0,1.2*Lmax,0.1*Lmax), minor=True)
ax.set_xticklabels([], minor=False)
ax.set_yticks(np.arange(ymin2,ymax2,dely), minor=False)
ax.set_yticklabels(np.arange(ymin2,ymax2,dely), minor=False, family='serif', size=14)
ax.set_yticks(np.arange(ymin2,ymax2,0.5*dely), minor=True)    
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.4f'))
plt.ylim(ymin2,ymax2)
plt.ylabel('$u_x$')
plt.xlim(0.0,Lmax)

for time in timearray:
  tt=float(time)
  filename='../output/xslice_t'+str(time)+'.dat'
  print('filename=',filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  rhoB=data[1]
  jB=data[2]
  epsilon=data[3]
  ux=data[4]
  pressure=data[5]
  ax.plot(x,ux,linestyle='-',linewidth=3,color='b')
  for i in range(0,len(x)):
    yguess[i]=0.001*sqrt(5.0/3.0)*sin(omega0*tt)*sin(k*x[i])
  plt.plot(x,yguess,linestyle='--',linewidth=3,color='k')

plt.savefig('xslice.pdf',format='pdf')
os.system('open -a Preview xslice.pdf')

quit()
