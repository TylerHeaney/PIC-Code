import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation
from simulator import Simulator
import random
import cProfile

matplotlib.use('QtAgg')


def two_stream():
  particle_num=40000
  cell_length=.125
  cell_number=400
  timestep=.06
  momentum=4e-26
  sd=momentum*.3


  half=(int)(particle_num/2)
  sim=Simulator(timestep,cell_length,cell_number,half*2,[1.67e-27]*half+[1.67e-27]*half,[1.6022e-19]*half+[1.6022e-19]*half,[1]*half*2)
  sim.initialize(0,0,[random.random()*cell_length*cell_number for _ in range(half*2)], [np.random.normal(momentum,sd) for _ in range(half)]+[np.random.normal(momentum,sd)*-1 for _ in range(half)])
  for i in range(len(sim.particles.momenta)):
    sim.particles.momenta[i] *= 1+.1*np.sin(2*np.pi*sim.particles.positions[i]/cell_number/cell_length)
  return sim


def rand():
  particle_num=40000
  cell_length=.125
  cell_number=400
  timestep=.06
  momentum=1e-26
  sd=0


  sim=Simulator(timestep,cell_length,cell_number,particle_num,[1.67e-27]*particle_num,[1.6022e-19]*particle_num,[1]*particle_num)
  sim.initialize(0,0,[random.random()*cell_length*cell_number for _ in range(particle_num)], [np.random.normal(momentum,sd) for _ in range(particle_num)])
  for i in range(len(sim.particles.momenta)):
    sim.particles.momenta[i] *= 1+.1*np.sin(2*np.pi*sim.particles.positions[i]/cell_number/cell_length)
  return sim

# Have same mass right now
def p_and_e():
  particle_num=40000
  cell_length=.125
  cell_number=400
  timestep=.06
  velocity=25
  sd=velocity*.3


  half=(int)(particle_num/2)
  sim=Simulator(timestep,cell_length,cell_number,half*2,[1.67e-27]*half+[1.67e-27]*half,[1.6022e-19]*half+[-1.6022e-19]*half,[1]*half*2)
  sim.initialize(0,0,[random.random()*cell_length*cell_number for _ in range(half*2)], [np.random.normal(velocity*1.67e-27,1.67e-27*sd) for _ in range(half)]+[np.random.normal(velocity*1.67e-27,1.67e-27*sd)*-1 for _ in range(half)])
  for i in range(len(sim.particles.momenta)):
    sim.particles.momenta[i] *= 1+.1*np.sin(2*np.pi*sim.particles.positions[i]/cell_number/cell_length)
  return sim

def p_and_e_in_both():
  particle_num=40000
  cell_length=.125
  cell_number=400
  timestep=.06
  velocity=25
  sd=velocity*.3


  half=(int)(particle_num/2)
  sim=Simulator(timestep,cell_length,cell_number,half*2,[1.67e-27]*half+[1.67e-27]*half,[1.6022e-19]*half+[-1.6022e-19]*half,[1]*half*2)
  sim.initialize(0,0,[random.random()*cell_length*cell_number for _ in range(half*2)], [np.random.normal(velocity*1.67e-27,1.67e-27*sd)*(-1)**_ for _ in range(half)]+[np.random.normal(velocity*1.67e-27,1.67e-27*sd)*(-1)**_ for _ in range(half)])
  for i in range(len(sim.particles.momenta)):
    sim.particles.momenta[i] *= 1+.1*np.sin(2*np.pi*sim.particles.positions[i]/cell_number/cell_length)
  return sim


def main():
  fig, ax = plt.subplots()
  ax.axhline(y=0, color='black', linewidth=0.5)
  choice = input("choose simulation:\n\t1) two-stream instability\n\t2) purely random\n\t3) protons and electrons\n\t4) p and e in both\n")
  sim = two_stream() if choice=='1' else rand() if choice=='2' else p_and_e() if choice=='3' else p_and_e_in_both()
  #cProfile.run('anim()')half=20
  plt.xticks([i*sim.delta_x for i in range(sim.num_cells)])

  update=anim(ax, sim)
  
  ani=FuncAnimation(fig,update,frames=500,interval=50)
  #ani.save('animation.gif', writer='imagemagick', fps=12)
  plt.show()




def anim(ax, sim):

  half=(int)(sim.particles.number / 2)

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
    p_p=sim.particles.velocities()[:half]
    x_e=sim.particles.positions[half:]
    p_e=sim.particles.velocities()[half:]
    ax.axhline(y=0, color='black', linewidth=0.5)
    plt.xticks([i*sim.delta_x for i in range(0,sim.num_cells,50)])
    protons,=ax.plot(x_p, p_p, 'o', markersize=.4, color='blue', label='Position: Protons')
    electrons,=ax.plot(x_e,p_e, 'o', markersize=.4, color='red', label='Position: Electrons')
    ax.set_ylim([-1e2,1e2])
    ax.set_xlim([0,sim.delta_x*sim.num_cells])

  return update



if(__name__=="__main__"):
  main()
  #cProfile.run('main()')