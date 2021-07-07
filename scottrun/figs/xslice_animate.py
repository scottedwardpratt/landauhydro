import matplotlib.animation as animation 
import matplotlib.pyplot as plt 
import numpy as np
import os

plt.figure(figsize=(6,10))
Lmax=100.0
delt=5
fig = plt.figure() 
nframes=1000
time=0
filename='../output/xslice_t'+str(time)+'.dat'
nframes=0
while os.path.isfile(filename):
  nframes+=1
  time=time+delt
  filename='../output/xslice_t'+str(time)+'.dat'

print('nframes=',nframes) 

axis = plt.axes(xlim =(0, Lmax), ylim=(0,1))
line, = axis.plot([],[],lw=2)
def init():
  ttl.set_text('')
  line.set_data([],[])
  return line,

ymin=0.0
ymax=0.6
axis1 = plt.axes(xlim =(0, Lmax), ylim=(ymin,ymax))
line1, = axis1.plot([],[],lw=2)
ttl=axis1.text(40,0.9*ymax+0.1*ymin,'',size='20')
def makedensity(i):
  time=i*delt
  tt=float(i*delt)
  filename='../output/xslice_t'+str(time)+'.dat'
  print('filename=',filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  #yguess=arange(0,len(x),1.0)
  rhoB=data[1]
  line1.set_data(x,rhoB)
  ttl.set_text('$\\rho_B$ vs. $x, t=$'+str(tt))
  return line1,
anim1=animation.FuncAnimation(fig, makedensity, init_func = init, frames = nframes, interval = 1, blit = True,save_count=None)
anim1.save('density.png',writer='Pillow',fps=15)
del anim1

ymin=0.0
ymax=0.3
axis2 = plt.axes(xlim =(0, Lmax),ylim=(ymin,ymax))
line2, = axis2.plot([],[],lw=2)
ttl=axis2.text(40,0.9*ymax+0.1*ymin,'',size='20')
def makeT(i):
  time=i*delt
  tt=float(i*delt)
  filename='../output/xslice_t'+str(time)+'.dat'
  print('filename=',filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  #yguess=arange(0,len(x),1.0)
  T=data[7]
  line2.set_data(x,T)
  ttl.set_text('$T$ vs. $x, t=$'+str(tt))
  return line2,
anim2=animation.FuncAnimation(fig, makeT, init_func = init, frames = nframes, interval = 1, blit = True,save_count=None)
anim2.save('T.png',writer='Pillow',fps=15)
del anim2

ymin=-0.025
ymax=0.025
axis3 = plt.axes(xlim =(0, Lmax),ylim=(ymin,ymax))
line3, = axis3.plot([],[],lw=2)
ttl=axis3.text(40,0.9*ymax+0.1*ymin,'',size='20')
def makeP(i):
  time=i*delt
  tt=float(i*delt)
  filename='../output/xslice_t'+str(time)+'.dat'
  print('filename=',filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  #yguess=arange(0,len(x),1.0)
  P=data[5]
  line3.set_data(x,P)
  ttl.set_text('$P$ vs. $x, t=$'+str(tt))
  return line3,
anim3=animation.FuncAnimation(fig, makeP, init_func = init, frames = nframes, interval = 1, blit = True,save_count=None)
anim3.save('P.png',writer='Pillow',fps=15)
del anim3

ymin=-0.015
ymax=0.05
axis4 = plt.axes(xlim =(0, Lmax),ylim=(ymin,ymax))
line4, = axis4.plot([],[],lw=2)
ttl=axis4.text(40,0.9*ymax+0.1*ymin,'',size='20')
def makeTxx(i):
  time=i*delt
  tt=float(i*delt)
  filename='../output/xslice_t'+str(time)+'.dat'
  print('filename=',filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  #yguess=arange(0,len(x),1.0)
  Txx=data[6]
  line4.set_data(x,Txx)
  ttl.set_text('$T_{xx}$ vs. $x, t=$'+str(tt))
  return line4,

anim4=animation.FuncAnimation(fig, makeTxx, init_func = init, frames = nframes, interval = 1, blit = True,save_count=None)
anim4.save('Txx.png',writer='Pillow',fps=15)
del anim4

quit()
