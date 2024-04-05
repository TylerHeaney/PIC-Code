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
  cell_length=.005
  cell_number=400
  timestep=.002
  momentum=5e-26
  sd=momentum*.1


  half=(int)(particle_num/2)
  sim=Simulator(timestep,cell_length,cell_number,half*2,[1.67e-27]*half+[1.67e-27]*half,[1.6022e-19]*half+[1.6022e-19]*half,[1]*half*2)
  sim.initialize(0,0,[random.random()*cell_length*cell_number for _ in range(half*2)], [np.random.normal(momentum,sd) for _ in range(half)]+[np.random.normal(0,sd)*-1 for _ in range(half)])
  for i in range(len(sim.particles.momenta)):
    sim.particles.momenta[i] *= 1+.1*np.sin(2*np.pi*sim.particles.positions[i]/cell_number/cell_length)
  return sim


def rand():
  particle_num=40000
  cell_length=.005
  cell_number=400
  timestep=.002
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


def oscil():
  particle_num=3
  cell_length=.00001
  cell_number=400
  timestep=.1
  velocity=0


  sim=Simulator(timestep,cell_length,cell_number,3,[1.67e-27]*1+[1.67]*2,[1.6022e-19]*1+[1.6022e-19]*2,[1]*3)
  sim.initialize(0,0,[cell_number*cell_length*.2,cell_length*cell_number*.1,cell_number*cell_length*.9], [0.0 for _ in range(1)]+[0.0 for _ in range(2)])
  return sim



def setup_plot(plt,fig,ax,sim):
  ax.cla()
  ax.axhline(y=0, color='black', linewidth=0.5)
  plt.xticks([i*sim.delta_x for i in range(0,sim.num_cells+1,50)])
  # ax.set_ylim([-1e2,1e2])
  ax.set_xlim([0,sim.delta_x*sim.num_cells])
  plt.xlabel("Position (cm)")
  plt.ylabel("Velocity (cm/sec)")

def setup_field_plot(plt,fig,ax,sim):
  plt.cla()
  plt.plot([i*sim.delta_x for i in range(len(sim.mesh.fields))],sim.mesh.fields,'-',markersize=.4)
  plt.xlabel("Position (cm)")
  plt.ylabel("Electric Field (V/cm)")
  # ax.set_ylim([-2e-5,2e-5])
  ax.set_xlim([0,sim.delta_x*sim.num_cells])



def plot_particles(plt,fig,ax,sim):
  half=(int)(sim.particles.number/2)
  x_p=sim.particles.positions[:half]
  p_p=sim.particles.velocities()[:half]
  x_e=sim.particles.positions[half:]
  p_e=sim.particles.velocities()[half:]
  # ax.plot(sim.particles.positions,sim.particles.fields,'o',markersize=.4)
  protons,=ax.plot(x_p, p_p, 'o', markersize=4, color='blue', label='Position: Protons')
  electrons,=ax.plot(x_e,p_e, 'o', markersize=4, color='red', label='Position: Electrons')

def plot_single_particle(plt,fig,ax,sim):
  global x,v
  jumps=[i for i in range(len(x)) if len(x)>0 and abs(x[i]-x[i-1])>1]
  jumps = jumps+[len(x)]
  # if len(jumps)>2:
    # jumps=jumps[-2:]
    # x=x[jumps[0]:]
  for n,j in enumerate(jumps):
    if n==0:
      ax.plot([i for i in x[:j]],[i for i in v[:j]],'-',color='black')
    else:
      ax.plot([i for i in x[jumps[n-1]:j]],[i for i in v[jumps[n-1]:j]],'-',color='black')


def main():
  fig, ax = plt.subplots()
  ax.axhline(y=0, color='black', linewidth=0.5)
  choice = input("choose simulation:\n\t1) two-stream instability\n\t2) purely random\n\t3) protons and electrons\n\t4) oscil\n")
  sim = two_stream() if choice=='1' else rand() if choice=='2' else p_and_e() if choice=='3' else oscil()
  #cProfile.run('anim()')half=20
  plt.xticks([i*sim.delta_x for i in range(sim.num_cells)])

  update=anim(fig,ax, sim)
  
  ani=FuncAnimation(fig,update,frames=500,interval=50)
  # ani.save('animation.gif', writer='imagemagick', fps=12)
  plt.show()
  fig, ax = plt.subplots()
  ax.cla()
  global x
  ax.plot(range(len(x)),[i[0] for i in x],color='blue')
  ax.plot(range(len(x)),[i[1] for i in x],color='red')
  plt.show()

  # x_p=sim.particles.positions[:half]
  # p_p=sim.particles.velocities()[:half]
  # x_e=sim.particles.positions[half:]
  # p_e=sim.particles.velocities()[half:]
  # ax.plot(sim.particles.positions,sim.particles.fields,'o',markersize=.4)
  # protons,=ax.plot(x_p, p_p, 'o', markersize=.4, color='blue', label='Position: Protons')
  # electrons,=ax.plot(x_e,p_e, 'o', markersize=.4, color='red', label='Position: Electrons')

def png():
  fig, ax = plt.subplots()
  ax.axhline(y=0, color='black', linewidth=0.5)
  sim = two_stream()
  plt.xticks([i*sim.delta_x for i in range(sim.num_cells)])

  half=(int)(sim.particles.number / 2)

  sim.step()
  
  global x

  setup_plot(plt,fig,ax,sim)
  plot_particles(plt,fig,ax,sim)
  # plot_single_particle(plt,fig,ax,sim)
  plt.savefig("t0.pdf",bbox_inches='tight')
  setup_field_plot(plt,fig,ax,sim)
  plt.savefig("f0.pdf",bbox_inches='tight')
  # plt.show()
  plt.cla()
  for i in range(25):
    sim.step()
    x.append([sim.particles.positions[0],sim.particles.velocities()[0]])
    # x.append([np.sum(.5*1.67e-27*np.abs(sim.particles.velocities()[half:])),np.sum(.5*1.67e-27*np.abs(sim.particles.velocities()[:half]))])

  
  setup_plot(plt,fig,ax,sim)
  plot_particles(plt,fig,ax,sim)
  # plot_single_particle(plt,fig,ax,sim)
  plt.savefig("t1.pdf",bbox_inches='tight')
  setup_field_plot(plt,fig,ax,sim)
  plt.savefig("f1.pdf",bbox_inches='tight')
  # plt.show()

  plt.cla()
  for i in range(25):
    sim.step()
    x.append([sim.particles.positions[0],sim.particles.velocities()[0]])
    # x.append([np.sum(.5*1.67e-27*np.abs(sim.particles.velocities()[half:])),np.sum(.5*1.67e-27*np.abs(sim.particles.velocities()[:half]))])
  
  setup_plot(plt,fig,ax,sim)

  plot_particles(plt,fig,ax,sim)
  # plot_single_particle(plt,fig,ax,sim)
  plt.savefig("t2.pdf",bbox_inches='tight')
  setup_field_plot(plt,fig,ax,sim)
  plt.savefig("f2.pdf",bbox_inches='tight')
  plt.show()



def anim(fig,ax, sim):

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
    
    global x,v,t
    # x.append([np.sum(.5*1.67e-27*np.abs(sim.particles.velocities()[half:])),np.sum(.5*1.67e-27*np.abs(sim.particles.velocities()[:half]))])
    x=np.append(x,[sim.particles.positions[0]])
    v=np.append(v,[sim.particles.velocities()[0]])
    t=np.append(t,sim.step_num*sim.delta_t)
    setup_plot(plt,fig,ax,sim)
    plot_particles(plt,fig,ax,sim)
    plot_single_particle(plt,fig,ax,sim)
    
    # ax.plot(sim.particles.positions,[0,0],'o',color='black')
    # ax.plot([x[i][0] for i in range(len(x)-10,len(x)) if i>0 and abs(x[i][0]-x[-1][0])<2 ],[x[i][1] for i in range(len(x)-10,len(x)) if i>0 and abs(x[i][0]-x[-1][0])<2 ],'-',markersize=.4)
    # if(len(x)==200):
      # plt.savefig("oscil.pdf")

  return update

x=np.array([])
v=np.array([])
t=np.array([])

if(__name__=="__main__"):
  main()
  # png()
  #cProfile.run('main()')