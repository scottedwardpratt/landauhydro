import matplotlib.animation as animation 
import matplotlib.pyplot as plt 
import numpy as np
import os

plt.figure(figsize=(6,10))
Lmax=100.0
delt=1
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
ymax=0.7
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

quit()
