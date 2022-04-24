import matplotlib.animation as animation 
import matplotlib.pyplot as plt 
import numpy as np
import os

plt.figure(figsize=(6,10))
Lmax=100.0
delt=5.0
fig1 = plt.figure() 
time=0
tstring="{:.1f}".format(time)
filename='../output/xslice_t'+tstring+'.dat'
nframes=0
while os.path.isfile(filename):
  nframes+=1
  time=time+delt
  tstring="{:.1f}".format(time)
  filename='../output/xslice_t'+tstring+'.dat'

print('nframes=',nframes) 

axis = plt.axes(xlim =(0, Lmax), ylim=(0,1))
line, = axis.plot([],[],lw=2)
def init1():
  ttl1.set_text('')
  line1.set_data([],[])
  return line1,

def init2():
  ttl2.set_text('')
  line2.set_data([],[])
  return line2,
  
def init3():
  ttl3.set_text('')
  line3.set_data([],[])
  return line3,

def init4():
  ttl4.set_text('')
  line4.set_data([],[])
  return line4,

ymin=0.0
ymax=0.7
axis1 = plt.axes(xlim =(0, Lmax), ylim=(ymin,ymax))
line1, = axis1.plot([],[],lw=2)
ttl1=axis1.text(40,0.9*ymax+0.1*ymin,'',size='20')
def makedensity(i):
  time=i*delt
  tt=float(i*delt)
  tstring="{:.1f}".format(time)
  filename='../output/xslice_t'+tstring+'.dat'
  print('filename=',filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  #yguess=arange(0,len(x),1.0)
  rhoB=data[1]
  line1.set_data(x,rhoB)
  ttl1.set_text('$\\rho_B$ vs. $x, t=$'+tstring)
  return line1,
anim1=animation.FuncAnimation(fig1, makedensity, init_func = init1, frames = nframes, interval = 1, blit = True,save_count=None)
anim1.save('density.png',writer='Pillow',fps=30)
del anim1

fig2 = plt.figure()
ymin=0.0
ymax=1.5
axis2 = plt.axes(xlim =(0, Lmax),ylim=(ymin,ymax))
line2, = axis2.plot([],[],lw=2)
ttl2=axis2.text(40,0.9*ymax+0.1*ymin,'',size='20')
def makeT(i):
  time=i*delt
  tt=float(i*delt)
  tstring="{:.1f}".format(time)
  filename='../output/xslice_t'+tstring+'.dat'
  print('filename=',filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  #yguess=arange(0,len(x),1.0)
  T=data[7]
  line2.set_data(x,T)
  ttl2.set_text('$T$ vs. $x, t=$'+tstring)
  return line2,
anim2=animation.FuncAnimation(fig2, makeT, init_func = init2, frames = nframes, interval = 1, blit = True,save_count=None)
anim2.save('T.png',writer='Pillow',fps=30)
del anim2

fig3 = plt.figure()
ymin=-1.0
ymax=1.0
axis3 = plt.axes(xlim =(0, Lmax),ylim=(ymin,ymax))
line3, = axis3.plot([],[],lw=2)
ttl3=axis3.text(40,0.9*ymax+0.1*ymin,'',size='20')
def makeP(i):
  time=i*delt
  tt=float(i*delt)
  tstring="{:.1f}".format(time)
  filename='../output/xslice_t'+tstring+'.dat'
  print('filename=',filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  #yguess=arange(0,len(x),1.0)
  P=data[5]
  line3.set_data(x,P)
  ttl3.set_text('$P$ vs. $x, t=$'+tstring)
  return line3,
anim3=animation.FuncAnimation(fig3, makeP, init_func = init3, frames = nframes, interval = 1, blit = True,save_count=None)
anim3.save('P.png',writer='Pillow',fps=30)
del anim3

fig4 = plt.figure()
ymin=-1.0
ymax=1.0
axis4 = plt.axes(xlim =(0, Lmax),ylim=(ymin,ymax))
line4, = axis4.plot([],[],lw=2)
ttl4=axis4.text(40,0.9*ymax+0.1*ymin,'',size='20')
def makeTxx(i):
  time=i*delt
  tt=float(i*delt)
  tstring="{:.1f}".format(time)
  filename='../output/xslice_t'+tstring+'.dat'
  print('filename=',filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  #yguess=arange(0,len(x),1.0)
  Txx=data[6]
  line4.set_data(x,Txx)
  ttl4.set_text('$T_{xx}$ vs. $x, t=$'+tstring)
  return line4,

anim4=animation.FuncAnimation(fig4, makeTxx, init_func = init4, frames = nframes, interval = 1, blit = True,save_count=None)
anim4.save('Txx.png',writer='Pillow',fps=30)
del anim4

quit()
