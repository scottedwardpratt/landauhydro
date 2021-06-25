import matplotlib.animation as animation 
import matplotlib.pyplot as plt 
import numpy as np
import os

plt.figure(figsize=(6,10))
Lmax=100.0
delt=2
fig = plt.figure() 
ymin=0.0
ymax=1.0
axis = plt.axes(xlim =(0, 100),
                ylim =(ymin,ymax) )
line, = axis.plot([],[],lw=2)
ttl=axis.text(80,0.9*ymax+0.1*ymin,'',size='20')
nframes=1000
time=0
filename='../output/xslice_t'+str(time)+'.dat'
nframes=0
while os.path.isfile(filename):
  nframes+=1
  time=time+delt
  filename='../output/xslice_t'+str(time)+'.dat'
  
print('nframes=',nframes) 

def init():
  ttl.set_text('')
  line.set_data([],[])
  return line,

delt=2
def makeline(i):
  time=i*delt
  tt=float(i*delt)
  filename='../output/xslice_t'+str(time)+'.dat'
  print('filename=',filename)
  data = np.loadtxt(filename,skiprows=0,unpack=True)
  x=data[0]
  #yguess=arange(0,len(x),1.0)
  rhoB=data[1]
  T=data[6]
  line.set_data(x,rhoB)
  ttl.set_text(str(tt))
  return line,

anim=animation.FuncAnimation(fig, makeline, init_func = init, frames = nframes, interval = 1, blit = True,save_count=None)
  
anim.save('density.png',writer='Pillow',fps=30)

quit()
