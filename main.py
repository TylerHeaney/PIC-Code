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
  cell_length=.025
  cell_number=400
  timestep=.06
  momentum=1e-26
  sd=momentum*.1


  half=(int)(particle_num/2)
  sim=Simulator(timestep,cell_length,cell_number,half*2,[1.67e-27]*half+[1.67e-27]*half,[1.6022e-19]*half+[1.6022e-19]*half,[1]*half*2)
  sim.initialize(0,0,[random.random()*cell_length*cell_number for _ in range(half*2)], [np.random.normal(momentum,sd) for _ in range(half)]+[np.random.normal(0,sd)*-1 for _ in range(half)])
  for i in range(len(sim.particles.momenta)):
    sim.particles.momenta[i] *= 1+.1*np.sin(2*np.pi*sim.particles.positions[i]/cell_number/cell_length)
  return sim


def rand():
  particle_num=40000
  cell_length=.025
  cell_number=400
  timestep=.06
  momentum=0
  sd=5e-27

  half=(int)(particle_num/2)
  sim=Simulator(timestep,cell_length,cell_number,particle_num,[1.67e-27]*particle_num,[1.6022e-19]*half+[-1.6022e-19]*half,[1]*particle_num)
  sim.initialize(0,0,[random.random()*cell_length*cell_number for _ in range(particle_num)], [np.random.normal(momentum,sd) for _ in range(particle_num)])
  # for i in range(len(sim.particles.momenta)):
    # sim.particles.momenta[i] *= 1+.1*np.sin(2*np.pi*sim.particles.positions[i]/cell_number/cell_length)
  return sim

# Have same mass right now
def p_and_e():
  particle_num=20000
  cell_length=.125
  cell_number=400
  timestep=.06
  velocity=5
  sd=velocity*.3


  half=(int)(particle_num/2)
  sim=Simulator(timestep,cell_length,cell_number,half*2,[1.67e-27]*half+[1.67e-31]*half,[1.6022e-19]*half+[-1.6022e-19]*half,[1]*half*2)
  sim.initialize(0,0,[random.random()*cell_length*cell_number for _ in range(half*2)], [np.random.normal(velocity*1.67e-27,1.67e-27*sd) for _ in range(half)]+[np.random.normal(velocity*1.67e-31,1.67e-31*sd)*-1 for _ in range(half)])
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
  # ani.save('animation.gif', writer='imagemagick', fps=12)
  plt.show()
  fig, ax = plt.subplots()
  ax.cla()
  global x
  ax.plot(range(len(x)),[i[0] for i in x],color='blue')
  ax.plot(range(len(x)),[i[1] for i in x],color='red')
  plt.show()


def png():
  fig, ax = plt.subplots()
  ax.axhline(y=0, color='black', linewidth=0.5)
  sim = two_stream()
  plt.xticks([i*sim.delta_x for i in range(sim.num_cells)])

  half=(int)(sim.particles.number / 2)

  sim.step()
  
  global x

  x_p=sim.particles.positions[:half]
  p_p=sim.particles.velocities()[:half]
  x_e=sim.particles.positions[half:]
  p_e=sim.particles.velocities()[half:]
  ax.axhline(y=0, color='black', linewidth=0.5)
  plt.xticks([i*sim.delta_x for i in range(0,sim.num_cells,50)])
  # ax.plot(sim.particles.positions,sim.particles.fields,'o',markersize=.4)
  protons,=ax.plot(x_p, p_p, 'o', markersize=.4, color='blue', label='Position: Protons')
  electrons,=ax.plot(x_e,p_e, 'o', markersize=.4, color='red', label='Position: Electrons')
  ax.set_ylim([-1e2,1e2])
  # ax.set_ylim([-2e-6,2e-6])
  ax.set_xlim([0,sim.delta_x*sim.num_cells])
  plt.savefig("t0.png",bbox_inches='tight')
  plt.cla()
  plt.plot([i*sim.delta_x for i in range(len(sim.mesh.fields))],sim.mesh.fields,'-',markersize=.4)
  ax.set_ylim([-2e-6,2e-6])
  plt.savefig("f0.png",bbox_inches='tight')
  # plt.show()
  plt.cla()
  for i in range(10):
    sim.step()
    x.append([np.sum(.5*1.67e-27*np.abs(sim.particles.velocities()[half:])),np.sum(.5*1.67e-27*np.abs(sim.particles.velocities()[:half]))])
  

  x_p=sim.particles.positions[:half]
  p_p=sim.particles.velocities()[:half]
  x_e=sim.particles.positions[half:]
  p_e=sim.particles.velocities()[half:]
  ax.axhline(y=0, color='black', linewidth=0.5)
  plt.xticks([i*sim.delta_x for i in range(0,sim.num_cells,50)])
  # ax.plot(sim.particles.positions,sim.particles.fields,'o',markersize=.4)
  protons,=ax.plot(x_p, p_p, 'o', markersize=.4, color='blue', label='Position: Protons')
  electrons,=ax.plot(x_e,p_e, 'o', markersize=.4, color='red', label='Position: Electrons')
  ax.set_ylim([-1e2,1e2])
  # ax.set_ylim([-2e-6,2e-6])
  ax.set_xlim([0,sim.delta_x*sim.num_cells])
  plt.savefig("t1.png",bbox_inches='tight')
  plt.cla()
  plt.plot([i*sim.delta_x for i in range(len(sim.mesh.fields))],sim.mesh.fields,'-',markersize=.4)
  ax.set_ylim([-2e-6,2e-6])
  plt.savefig("f1.png",bbox_inches='tight')
  # plt.show()

  plt.cla()
  for i in range(90):
    sim.step()
    x.append([np.sum(.5*1.67e-27*np.abs(sim.particles.velocities()[half:])),np.sum(.5*1.67e-27*np.abs(sim.particles.velocities()[:half]))])
  
  x_p=sim.particles.positions[:half]
  p_p=sim.particles.velocities()[:half]
  x_e=sim.particles.positions[half:]
  p_e=sim.particles.velocities()[half:]
  ax.axhline(y=0, color='black', linewidth=0.5)
  plt.xticks([i*sim.delta_x for i in range(0,sim.num_cells,50)])
  # ax.plot(sim.particles.positions,sim.particles.fields,'o',markersize=.4)
  protons,=ax.plot(x_p, p_p, 'o', markersize=.4, color='blue', label='Position: Protons')
  electrons,=ax.plot(x_e,p_e, 'o', markersize=.4, color='red', label='Position: Electrons')
  ax.set_ylim([-1e2,1e2])
  # ax.set_ylim([-2e-6,2e-6])
  ax.set_xlim([0,sim.delta_x*sim.num_cells])
  plt.savefig("t2.png",bbox_inches='tight')
  plt.cla()
  plt.plot([i*sim.delta_x for i in range(len(sim.mesh.fields))],sim.mesh.fields,'-',markersize=.4)
  ax.set_ylim([-2e-6,2e-6])
  plt.savefig("f2.png",bbox_inches='tight')
  plt.show()



def anim(ax, sim):

  half=(int)(sim.particles.number / 2)

  x_p=sim.particles.positions[:half]
  p_p=sim.particles.velocities()[:half]
  x_e=sim.particles.positions[half:]
  p_e=sim.particles.velocities()[half:]
  protons,=ax.plot(x_p, p_p, 'o', markersize=.4, color='blue', label='Position: Protons')
  electrons,=ax.plot(x_e,p_e, 'o', markersize=.4, color='red', label='Position: Electrons')

  

  def update(frame):
    for i in range(1):
      sim.step()
    ax.cla()
    
    global x
    # x.append([np.sum(.5*1.67e-27*np.abs(sim.particles.velocities()[half:])),np.sum(.5*1.67e-27*np.abs(sim.particles.velocities()[:half]))])
    x.append([sim.particles.positions[0],sim.particles.velocities()[0]])

    x_p=sim.particles.positions[:half]
    p_p=sim.particles.velocities()[:half]
    x_e=sim.particles.positions[half:]
    p_e=sim.particles.velocities()[half:]
    ax.axhline(y=0, color='black', linewidth=0.5)
    plt.xticks([i*sim.delta_x for i in range(0,sim.num_cells,50)])
    # ax.plot(sim.particles.positions,sim.particles.fields,'o',markersize=.4)
    protons,=ax.plot(x_p, p_p, 'o', markersize=.4, color='blue', label='Position: Protons')
    electrons,=ax.plot(x_e,p_e, 'o', markersize=.4, color='red', label='Position: Electrons')
    ax.plot([i[0] for i in x],[i[1] for i in x],'-',markersize=.4)
    ax.set_ylim([-1e1,1e1])
    # ax.set_ylim([-2e-6,2e-6])
    ax.set_xlim([0,sim.delta_x*sim.num_cells])

  return update

x=[]

if(__name__=="__main__"):
  main()
  # png()
  #cProfile.run('main()')