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

modeldir='../../output_posterior_corrected/default/'

ipanel=0
ax0 = fig.add_axes([xx0,yy0+ipanel*hh0/Npanels,ww0,hh0/Npanels])
ipanel=1
ax1 = fig.add_axes([xx0,yy0+ipanel*hh0/Npanels,ww0,hh0/Npanels])
ipanel=2
ax2 = fig.add_axes([xx0,yy0+ipanel*hh0/Npanels,ww0,hh0/Npanels])
ax0.set_xticks(np.arange(0,1.2*Lmax,0.2*Lmax), minor=False)
ax0.set_xticks(np.arange(0,1.2*Lmax,0.1*Lmax), minor=True)
ax0.set_xticklabels(np.arange(0,1.2*Lmax,0.2*Lmax), minor=False, size=14)
ax0.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
#ax0.xlabel('$x$',fontsize=24)
ax1.set_xticks(np.arange(0,1.2*Lmax,0.2*Lmax), minor=False)
ax1.set_xticks(np.arange(0,1.2*Lmax,0.1*Lmax), minor=True)
ax1.set_xticklabels([], minor=False)
ax2.set_xticks(np.arange(0,1.2*Lmax,0.2*Lmax), minor=False)
ax2.set_xticks(np.arange(0,1.2*Lmax,0.1*Lmax), minor=True)
ax2.set_xticklabels([], minor=False)

# ipanel=0
ymax0=1.015
ymin0=0.985
dely=0.005
ax0.set_yticks(np.arange(ymin0,1.1*ymax0,dely), minor=False)
ax0.set_yticklabels(np.arange(ymin0,1.1*ymax0,dely), minor=False, family='serif', size=14)
ax0.set_yticks(np.arange(ymin0,1.1*ymax0,0.5*dely), minor=True)
ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
#ax0.ylabel('$\\rho$')

# ipanel=1
ymax1=0.97
ymin1=1.03
dely=0.01
ax1.set_yticks(np.arange(ymin1,1.1*ymax1,dely), minor=False)
ax1.set_yticklabels(np.arange(ymin1,1.1*ymax1,dely), minor=False, family='serif', size=14)
ax1.set_yticks(np.arange(ymin1,1.1*ymax1,0.5*dely), minor=True)
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
#ax1.ylabel('$P$')

# ipanel=2
ymax2=0.02
ymin2=-ymax2
dely=0.005
ax2.set_yticks(np.arange(ymin2,1.1*ymax2,dely), minor=False)
ax2.set_yticklabels(np.arange(ymin2,1.1*ymax2,dely), minor=False, family='serif', size=14)
ax2.set_yticks(np.arange(ymin2,1.1*ymax2,0.5*dely), minor=True)    
ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.4f'))

#for time in ['200','220','240','260','280','300']:
for time in ['200']:
  tt=float(time)
  print('tt=',tt)
  filename='../output/xslice_t'+str(time)+'.dat'
  print('filename=',filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  yguess=x
  rhoB=data[1]
  jB=data[2]
  epsilon=data[3]
  ux=data[4]
  pressure=data[5]
  
  #plt.ylabel("\\rho")
  ax0.plot(x,rhoB,linestyle='-',linewidth=3,color='r')
      #plt.clf()
  for i in range(0,len(x)):
    yguess[i]=1.0+0.01*cos(omega0*tt)*cos(k*x[i])
  #ax0.plot(x,yguess,linestyle='--',linewidth=3,color='k')
  
  #plt.ylim(ymin1,ymax1)
  #plt.ylabel('$P')
  ax1.plot(x,pressure,linestyle='-',linewidth=3,color='g')
  
  #plt.ylim(ymin2,ymax2)
  #plt.ylabel('$u_x$')
  ax2.plot(x,ux,linestyle='-',linewidth=3,color='b')

  #plt.xlim(0.0,Lmax)
  #plt.ylim(ymin0,ymax0)

    
 



#plt.savefig('xslice.pdf',format='pdf')
#os.system('open -a Preview xslice.pdf')
plt.show()

quit()
