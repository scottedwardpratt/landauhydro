import matplotlib.animation as animation 
import matplotlib.pyplot as plt 
import numpy as np
import os

Idata=[1,5,6,7]
Lmax=100.0
delt=5.0
data_labels=['$x$','$\\rho_B$','$J_x$','$\\epsilon$','$P$','$Txx$','$T$','$K_x$','$c_s^2$']
ifig=int(input('input 0 for rho, 1 for T, 2 for P or 3 for Txx: '))

plt.figure(figsize=(3,5))

fig = plt.figure() 
time=0
tstring="{:06.1f}".format(time)
filename='../output/xslice_t'+tstring+'.dat'
#print('reading from '+filename)
nframes=0
while os.path.isfile(filename):
  nframes+=1
  time=time+delt
  tstring="{:06.1f}".format(time)
  filename='../output/xslice_t'+tstring+'.dat'
  print('reading from '+filename)
print('nframes=',nframes) 

ymin=0.0
ymax=0.7
axis = plt.axes(xlim =(0, Lmax), ylim=(ymin,ymax))
ttl=axis.text(40,0.9*ymax+0.1*ymin,'',size='20')
line, = axis.plot([],[],lw=2)
def init():
  ttl.set_text('')
  line.set_data([],[])
  return line,

def makedensity(i):
  time=i*delt
  tt=float(i*delt)
  tstring="{:06.1f}".format(time)
  filename='../output/xslice_t'+tstring+'.dat'
  print('reading '+filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  #yguess=arange(0,len(x),1.0)
  rhoB=data[1]
  line.set_data(x,rhoB)
  ttl.set_text('$\\rho_B$ vs. $x, t=$'+tstring)
  return line,
def makeT(i):
  time=i*delt
  tt=float(i*delt)
  tstring="{:06.1f}".format(time)
  filename='../output/xslice_t'+tstring+'.dat'
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  #yguess=arange(0,len(x),1.0)
  T=data[7]
  line.set_data(x,T)
  ttl.set_text('$T$ vs. $x, t=$'+tstring)
  return line,
def makeP(i):
  time=i*delt
  tt=float(i*delt)
  tstring="{:06.1f}".format(time)
  filename='../output/xslice_t'+tstring+'.dat'
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  #yguess=arange(0,len(x),1.0)
  P=data[5]
  line.set_data(x,P)
  ttl.set_text('$P$ vs. $x, t=$'+tstring)
  return line,
def makeTxx(i):
  time=i*delt
  tt=float(i*delt)
  tstring="{:06.1f}".format(time)
  filename='../output/xslice_t'+tstring+'.dat'
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  #yguess=arange(0,len(x),1.0)
  Txx=data[6]
  line.set_data(x,Txx)
  ttl.set_text('$T_{xx}$ vs. $x, t=$'+tstring)
  delete(data)
  return line,


if ifig==0:
  anim=animation.FuncAnimation(fig, makedensity, init_func = init, frames = nframes, interval = 1, blit = True,save_count=0)
  anim.save('density.png',writer='Pillow',fps=30)
  del anim
if ifig==1:
  anim=animation.FuncAnimation(fig, makeT, init_func = init, frames = nframes, interval = 1, blit = True,save_count=0)
  anim.save('T.png',writer='Pillow',fps=30)
  del anim
if ifig==2:
  anim=animation.FuncAnimation(fig, makeP, init_func = init, frames = nframes, interval = 1, blit = True,save_count=0)
  anim.save('T.png',writer='Pillow',fps=30)
  del anim
if ifig==3:
  anim=animation.FuncAnimation(fig, makeTxx, init_func = init, frames = nframes, interval = 1, blit = True,save_count=0)
  anim.save('density.png',writer='Pillow',fps=30)
  del anim


quit()
