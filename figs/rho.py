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

#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
plt.rc('text', usetex=False)
plt.figure(figsize=(6,5))
fig = plt.figure(1)
ax = fig.add_axes([0.15,0.12,0.8,0.8])

for t in range(10,100,10):
  filename='../run/output/current_t'+str(t)+'.dat'
  print 'filename=',filename
  mydata = np.loadtxt(filename,skiprows=0,unpack=True)

  x=mydata[0]
  rho=mydata[5]
  plt.plot(x,rho,linestyle='-',linewidth=2,color='r')


ax.tick_params(axis='both', which='major', labelsize=14)

ax.set_xticks(np.arange(0,120,20), minor=False)
ax.set_xticklabels(np.arange(120,41,20), minor=False, family='serif')
ax.set_xticks(np.arange(0,121,10), minor=True)
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
plt.xlim(0.0,100)

ax.set_yticks(np.arange(0,3,0.5), minor=False)
ax.set_yticklabels(np.arange(0,3,0.5), minor=False, family='serif')
ax.set_yticks(np.arange(0,3,0.1), minor=True)
plt.ylim(0.0,2.0)
#ax.set_yticks(0.1:1.0:10.0:100.0, minor=True)
#ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
ax.yaxis.set_major_formatter(sformatter)

plt.xlabel('$t$ (s)', fontsize=18, weight='normal')
plt.ylabel('voltage=$RI$ (mV)',fontsize=18)
#plt.title('MathText Number $\sum_{n=1}^\infty({-e^{i\pi}}/{2^n})$!',
#fontsize=12, color='gray')
#plt.subplots_adjust(top=0.85)
plt.savefig('rho.pdf',format='pdf')
os.system('open -a Preview rho.pdf')
quit()
