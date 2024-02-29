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

    fig, ax = plt.subplots()
    ax.scatter(positions, [0]*len(positions), color='blue', label='Position')
    
    for i, particle in enumerate(particles):
      print(f'{i}: pos: {particle.position} | momentum: {particle.momentum}')
      
      
    ax.axhline(y=0, color='black', linewidth=0.5)
    plt.xticks([i*delta_x for i in range(num_cells)])
    ax.set_xlabel('Position')
    ax.set_title('Particles on a Number Line')
    ax.legend()
    plt.show()


def main():
   sim=Simulator(1,.1,5,5,1,1,1)
   sim.initialize(0,0,[(.1,.1),(.26,0),(.44,0),(.3,0),(.15,0)])

   plot(sim.mesh, sim.particles,sim.delta_x,sim.num_cells)
   sim.step()
   sim.step()
   plot(sim.mesh, sim.particles,sim.delta_x,sim.num_cells)


if(__name__=="__main__"):
    main()