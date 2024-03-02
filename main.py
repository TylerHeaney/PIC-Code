import numpy
import matplotlib.pyplot as plt
import matplotlib
import particles
import mesh
from simulator import Simulator

matplotlib.use('QtAgg')

def plot(mesh, particles,delta_x,num_cells):
  positions = [particle.position for particle in particles]
  momenta = [particle.momentum for particle in particles]
  field=[n.field for n in mesh]
  fig, ax = plt.subplots()
  ax.scatter(positions[:20], [0]*(int)(len(positions)/2), color='blue', label='Position: Protons')
  ax.scatter(positions[20:],[0]*(int)(len(positions)/2), color='red', label='Position: Electrons')
  ax.plot(numpy.linspace(0,10,len(mesh)),field,linewidth=2)
  for i, particle in enumerate(particles):
    print(f'{i}: pos: {particle.position} | momentum: {particle.momentum}')
      
      
  ax.axhline(y=0, color='black', linewidth=0.5)
  plt.xticks([i*delta_x for i in range(num_cells)])
  ax.set_xlabel('Position')
  ax.legend()
  plt.show()


def main():
  sim=Simulator(.01,.1,100,40,[1.67e-27]*20+[9.11e-31]*20,[1.6022e-19]*20+[-1.6022e-19]*20,[1]*40)
  sim.initialize(0,0,[(.9,0),(1.0,0),(1.1,0)])
   
  while(1):
    plot(sim.mesh, sim.particles,sim.delta_x,sim.num_cells)
    sim.step()
    #print([n.field for n in sim.mesh])
    print([p.momentum for p in sim.particles])


if(__name__=="__main__"):
  main()