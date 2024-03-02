import numpy
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation
import particles
import mesh
from simulator import Simulator
import random

matplotlib.use('QtAgg')

def plot(mesh, particles,delta_x,num_cells):
  positions = [particle.position for particle in particles]
  momenta = [particle.momentum for particle in particles]
  fig, ax = plt.subplots()
  
  for i, particle in enumerate(particles):
    print(f'{i}: pos: {particle.position} | momentum: {particle.momentum}')
      
      
  ax.set_xlabel('Position')
  ax.legend()
  plt.show()

def anim():
  half=20
  sim=Simulator(.001,.1,100,40,[1.67e-27]*half+[9.11e-31]*half,[1.6022e-19]*half+[-1.6022e-19]*half,[1]*half*2)
  sim.initialize(0,0,[(.9,0),(1.0,0),(1.1,0)])

  fig, ax = plt.subplots()
  ax.axhline(y=0, color='black', linewidth=0.5)
  plt.xticks([i*sim.delta_x for i in range(sim.num_cells)])

  x_p=[p.position for p in sim.particles[:half]]
  x_e=[e.position for e in sim.particles[half:]]
  protons,=ax.plot(x_p, [0]*(len(x_p)), 'o', color='blue', label='Position: Protons')
  electrons,=ax.plot(x_e,[0]*(len(x_e)), 'o', color='red', label='Position: Electrons')
  field=[n.field for n in sim.mesh]
  f,=ax.plot(numpy.linspace(0,10,len(sim.mesh)),field,linewidth=2)

  def update(frame):
    for i in range(2):
      sim.step()
    ax.cla()
    
    x_p=[p.position for p in sim.particles[:half]]
    x_e=[e.position for e in sim.particles[half:]]
    ax.axhline(y=0, color='black', linewidth=0.5)
    plt.xticks([i*sim.delta_x for i in range(sim.num_cells)])
    ax.set_ylim([-1e-8,1e-8])
    protons,=ax.plot(x_p, [0]*(len(x_p)), 'o', color='blue', label='Position: Protons')
    electrons,=ax.plot(x_e,[0]*(len(x_e)), 'o', color='red', label='Position: Electrons')
    field=[n.field for n in sim.mesh]
    f,=ax.plot(numpy.linspace(0,10,len(sim.mesh)),field,linewidth=2)
    # protons.set_xdata([p.position for p in sim.particles[:half]])
    # electrons.set_xdata([e.position for e in sim.particles[half:]])
    # f.set_ydata([n.field for n in sim.mesh])
    #return protons,electrons,f,
  
  anim=FuncAnimation(fig,update,frames=500,interval=50)
  anim.save('animation.gif', writer='imagemagick', fps=12)
  #plt.show()

def main():
  anim()
  # while(1):
  #   sim.step()
  #   #print([n.field for n in sim.mesh])
  #   print([p.momentum for p in sim.particles])
  #   anim(sim.particles)


if(__name__=="__main__"):
  main()