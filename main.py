import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation
from simulator import Simulator
import random
import cProfile

matplotlib.use('QtAgg')

def anim():
  particle_num=20000
  cell_length=.125
  cell_number=400
  timestep=.06
  momentum=4e-26
  sd=momentum*.3


  half=(int)(particle_num/2)
  sim=Simulator(timestep,cell_length,cell_number,half*2,[1.67e-27]*half+[1.67e-27]*half,[1.6022e-19]*half+[1.6022e-19]*half,[1]*half*2)
  sim.initialize(0,0,[random.random()*cell_length*cell_number for _ in range(half*2)], [np.random.normal(momentum,sd) for _ in range(half)]+[np.random.normal(momentum,sd)*-1 for _ in range(half)])
  fig, ax = plt.subplots()
  for i in range(len(sim.particles.momenta)):
    sim.particles.momenta[i] *= 1+.1*np.sin(2*np.pi*sim.particles.positions[i]/cell_number)
  ax.axhline(y=0, color='black', linewidth=0.5)
  plt.xticks([i*sim.delta_x for i in range(sim.num_cells)])

  x_p=sim.particles.positions[:half]
  p_p=sim.particles.momenta[:half]
  x_e=sim.particles.positions[half:]
  p_e=sim.particles.momenta[half:]
  protons,=ax.plot(x_p, p_p, 'o', markersize=.4, color='blue', label='Position: Protons')
  electrons,=ax.plot(x_e,p_e, 'o', markersize=.4, color='red', label='Position: Electrons')

  def update(frame):
    for i in range(1):
      sim.step()
    ax.cla()
    
    x_p=sim.particles.positions[:half]
    p_p=sim.particles.momenta[:half]
    x_e=sim.particles.positions[half:]
    p_e=sim.particles.momenta[half:]
    ax.axhline(y=0, color='black', linewidth=0.5)
    plt.xticks([i*sim.delta_x for i in range(0,sim.num_cells,50)])
    protons,=ax.plot(x_p, p_p, 'o', markersize=.4, color='blue', label='Position: Protons')
    electrons,=ax.plot(x_e,p_e, 'o', markersize=.4, color='red', label='Position: Electrons')
    ax.set_ylim([-1e-25,1e-25])
    ax.set_xlim([0,cell_length*cell_number])

  anim=FuncAnimation(fig,update,frames=100,interval=50)
  anim.save('animation.gif', writer='imagemagick', fps=12)

def main():
  anim()
  #cProfile.run('anim()')half=20

if(__name__=="__main__"):
  main()
  #cProfile.run('main()')