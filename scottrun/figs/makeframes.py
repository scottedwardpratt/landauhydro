import matplotlib.animation as animation 
import matplotlib.pyplot as plt 
import numpy as np
import os

from matplotlib import ticker
from matplotlib.ticker import ScalarFormatter
#sformatter=ScalarFormatter(useOffset=True,useMathText=True)
#sformatter.set_scientific(True)
#sformatter.set_powerlimits((-2,3))

delt=5
xmax=100
nframes=0
filename='../output/xslice_t'+"0000.0"+'.dat'
print('reading from '+filename)
time=0.0
while os.path.isfile(filename):
  nframes+=1
  time=time+delt
  tstring="{:06.1f}".format(time)
  filename='../output/xslice_t'+tstring+'.dat'
print('nframes=',nframes)

plt.figure(figsize=(5,9))
fig=plt.figure(1)
bottom=0.08
height=0.25*(0.99-bottom)
x0=0.175
width=0.94-x0

for iframe in range(0,nframes+2):
  time=iframe*delt
  tstring="{:06.1f}".format(time)
  filename='../output/xslice_t'+tstring+'.dat'
  print('reading '+filename)
  mydata=np.loadtxt(filename,skiprows=0,unpack=True)
  #x=mydata[0]
  #yguess=arange(0,len(x),1.0)
  #rhoB=mydata[1]
  #T=mydata[7]
  #P=mydata[5]
  #Txx=mydata[6]
 
  ################## Panel 1 #####################

  ymin=-2.0+0.0001
  ymax=2.0
  ax = fig.add_axes([x0,bottom,width,height])
  plt.plot(mydata[0],mydata[1],linestyle='-',color='r',markersize=6,marker=None)
  ax.set_xticks(np.arange(0,xmax+0.001,xmax/5), minor=False)
  ax.set_xticklabels(np.arange(0,xmax+0.001,xmax/5), minor=False, family='sans',fontsize=18)
  ax.set_xticks(np.arange(0,xmax+0.001,xmax/10), minor=True)
  ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
  ax.set_yticks(np.arange(ymin,ymax,0.2), minor=False)
  ax.set_yticklabels(np.arange(ymin,ymax,0.2), minor=False, family='sans',fontsize=18)
  ax.set_yticks(np.arange(ymin,ymax,0.1), minor=True)
  ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
  plt.xlim(0,xmax)
  plt.ylim(0,0.7)
  plt.ylabel("$\\rho$",fontsize=24,labelpad=-1)
  plt.xlabel("$x$",fontsize=24)

  ################## Panel 2 #####################

  ymin=-2.0+0.0001
  ymax=2.0
  ax = fig.add_axes([x0,bottom+height,width,height])
  plt.plot(mydata[0],mydata[7],linestyle='-',color='r',markersize=6,marker=None)
  ax.set_xticks(np.arange(0,xmax,xmax/5), minor=False)
  ax.set_xticklabels([], minor=False, family='sans',fontsize=18)
  ax.set_xticks(np.arange(0,xmax,xmax/10), minor=True)
  ax.set_yticks(np.arange(ymin,ymax,0.2), minor=False)
  ax.set_yticklabels(np.arange(ymin,ymax,0.2), minor=False, family='sans',fontsize=18)
  ax.set_yticks(np.arange(ymin,ymax,0.1), minor=True)
  ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
  plt.xlim(0.0,xmax)
  plt.ylim(0.0,0.3)
  plt.ylabel("$T$",fontsize=24,labelpad=-3)
  plt.xlabel(None)

  ################## Panel 3 #####################

  ymin=-2.0+0.0001
  ymax=2.0
  ax = fig.add_axes([x0,bottom+2*height,width,height])
  plt.plot(mydata[0],mydata[5],linestyle='-',color='r',markersize=6,marker=None)
  plt.xlim(0.0,100)
  plt.ylim(-0.2,0.4)
  ax.set_xticks(np.arange(0,xmax,xmax/5), minor=False)
  ax.set_xticklabels([], minor=False, family='sans',fontsize=18)
  ax.set_xticks(np.arange(0,xmax,xmax/10), minor=True)
  ax.set_yticks(np.arange(ymin,ymax,0.2), minor=False)
  ax.set_yticklabels(np.arange(ymin,ymax,0.2), minor=False, family='sans',fontsize=18)
  ax.set_yticks(np.arange(ymin,ymax,0.1), minor=True)
  ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
  plt.xlim(0.0,xmax)
  plt.ylim(-0.15,0.3)
  plt.ylabel("$P$",fontsize=24,labelpad=-3)
  plt.xlabel(None)

################## Panel 4 #####################

  ymin=-2.0+0.0001
  ymax=2.0
  ax = fig.add_axes([x0,bottom+3*height,width,height])
  plt.plot(mydata[0],mydata[6],linestyle='-',color='r',markersize=6,marker=None)
  ax.set_xticks(np.arange(0,xmax,xmax/5), minor=False)
  ax.set_xticklabels([], minor=False, family='sans',fontsize=18)
  ax.set_xticks(np.arange(0,xmax,xmax/10), minor=True)
  ax.set_yticks(np.arange(ymin,ymax,0.2), minor=False)
  ax.set_yticklabels(np.arange(ymin,ymax,0.2), minor=False, family='sans',fontsize=18)
  ax.set_yticks(np.arange(ymin,ymax,0.1), minor=True)
  ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
  plt.xlim(0.0,xmax)
  plt.ylim(-0.15,0.3)
  plt.ylabel("$T_{xx}$",fontsize=24,labelpad=-3)
  plt.xlabel(None)


  plt.text(60,0.21,"$t=$"+str(time),fontsize=24,ha='left')

###############################################

  tstring="{:06.1f}".format(time)
  filename="pdffiles/t"+tstring+".pdf"
  plt.savefig(filename)
  plt.cla()
  plt.clf()
  del mydata


quit()
