import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation
import particles
import mesh
from simulator import Simulator
import random
import cProfile

matplotlib.use('QtAgg')

def anim():
  particle_num=500
  cell_length=.125
  cell_number=400
  timestep=.125
  momentum=1e-26
  sd=momentum*.3



  half=(int)(particle_num/2)
  sim=Simulator(timestep,cell_length,cell_number,half*2,[1.67e-27]*half+[1.67e-27]*half,[1.6022e-19]*half+[-1.6022e-19]*half,[1]*half*2)
  # sim.initialize(0,0,[])
  sim.initialize(0,0,[(random.random()*cell_length*cell_number,np.random.normal(momentum,sd)) for _ in range(half)]+[(random.random()*cell_length*cell_number,-1*np.random.normal(momentum,sd)) for _ in range(half)])
  fig, ax = plt.subplots()
  for p in sim.particles:
    p.momentum *= 1+.1*np.sin(2*np.pi*p.position/cell_number)
  ax.axhline(y=0, color='black', linewidth=0.5)
  plt.xticks([i*sim.delta_x for i in range(sim.num_cells)])

  x_p=[p.position for p in sim.particles[:half]]
  p_p=[p.momentum for p in sim.particles[:half]]
  x_e=[e.position for e in sim.particles[half:]]
  p_e=[e.momentum for e in sim.particles[half:]]
  print(p_p)
  protons,=ax.plot(x_p, p_p, 'o', color='blue', label='Position: Protons')
  electrons,=ax.plot(x_e,p_e, 'o', color='red', label='Position: Electrons')
  #field=[n.field for n in sim.mesh]
  #f,=ax.plot(numpy.linspace(0,10,len(sim.mesh)),field,linewidth=2)

  def update(frame):
    for i in range(2):
      sim.step()
    ax.cla()
    
    x_p=[p.position for p in sim.particles[:half]]
    p_p=[p.momentum for p in sim.particles[:half]]
    x_e=[e.position for e in sim.particles[half:]]
    p_e=[e.momentum for e in sim.particles[half:]]
    ax.axhline(y=0, color='black', linewidth=0.5)
    plt.xticks([i*sim.delta_x for i in range(0,sim.num_cells,50)])
    protons,=ax.plot(x_p, p_p, 'o', color='blue', label='Position: Protons')
    electrons,=ax.plot(x_e,p_e, 'o', color='red', label='Position: Electrons')
    ax.set_ylim([-1e-25,1e-25])
    ax.set_xlim([0,cell_length*cell_number])
    #protons,=ax.plot(x_p, [0]*(len(x_p)), 'o', color='blue', label='Position: Protons')
    #electrons,=ax.plot(x_e,[0]*(len(x_e)), 'o', color='red', label='Position: Electrons')
    #field=[n.field for n in sim.mesh]
    #f,=ax.plot(numpy.linspace(0,10,len(sim.mesh)),field,linewidth=2)
  
  anim=FuncAnimation(fig,update,frames=50,interval=50)
  anim.save('animation.gif', writer='imagemagick', fps=12)
  #plt.show()

def main():
  anim()
  #cProfile.run('anim()')half=20
  # half=20
  # sim=Simulator(.001,.1,100,half*2,[1.67e-27]*half+[9.11e-31]*half,[1.6022e-19]*half+[-1.6022e-19]*half,[1]*half*2)
  # sim.initialize(0,0,[(.9,0),(1.0,0),(1.1,0)])

  # fig, ax = plt.subplots()
  # ax.axhline(y=0, color='black', linewidth=0.5)
  # plt.xticks([i*sim.delta_x for i in range(sim.num_cells)])

  # x_p=[p.position for p in sim.particles[:half]]
  # x_e=[e.position for e in sim.particles[half:]]
  # protons,=ax.plot(x_p, [0]*(len(x_p)), 'o', color='blue', label='Position: Protons')
  # electrons,=ax.plot(x_e,[0]*(len(x_e)), 'o', color='red', label='Position: Electrons')
  # field=[n.field for n in sim.mesh]
  # f,=ax.plot(numpy.linspace(0,10,len(sim.mesh)),field,linewidth=2)
  # for i in range(100):
  #   sim.step()


if(__name__=="__main__"):
  main()
  #cProfile.run('main()')