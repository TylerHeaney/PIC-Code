import numpy as np
from particles import Particle
from mesh import Node
import random

class Simulator:
    """The actual simulator"""

    def __init__(self,delta_t,delta_x,num_cells,particle_num,particle_mass,particle_charge,particle_agg):
        self.step_num=0
        self.delta_t=delta_t
        self.delta_x=delta_x
        self.num_cells=num_cells
        self.create_particles(particle_num,particle_agg,particle_mass,particle_charge)
        self.create_mesh(delta_x,num_cells)

    def initialize(self, lhw, rhw, positions_and_momentums):
        if(len(positions_and_momentums)!=len(self.particles)): # use random numbers
            for p in self.particles:
                p.initialize(random.random()*self.delta_x*self.num_cells,random.random())
        else:
            for i in range(len(self.particles)):
                self.particles[i].initialize(*positions_and_momentums[i])
        
        self.mesh[0].charge=lhw
        self.mesh[-1].charge=rhw

    def step(self):
        self.scatter()
        self.field_solve()
        self.gather()
        self.push()


        self.step_num+=1
   
    def zero_charge(self):
        for node in self.mesh:
            node.charge=0
   

    def scatter(self):
        self.zero_charge()
        for node in self.mesh:
            node.aggregate_charge(self.particles)
    
    def zero_field(self):
        for p in self.particles:
            p.field=0
            
    def mesh_zero_field(self):
        for node in mesh:
            node.field=0
            

    def field_solve(self):
        self.mesh_zero_field()
        A=self.create_poisson_matrix
        potential_vector=np.linalg.solve(A,self.mesh.charge_vector)
        self.mesh.find_E(potential_vector)
    
    def gather(self):
        self.zero_field()
        for particle in self.particles:
            particle.aggregate_field(self.mesh)
    

    def time_step(self):
        return self.step_num*self.delta_t

    def push(self):
        for particle in self.particles:
            particle.position=particle.position+particle.velocity()*self.time_step()
            particle.momentum=particle.momentum+particle.charge*particle.field*self.time_step()


        
    def create_particles(self,particle_num,particle_agg,particle_mass,particle_charge):
        self.particles=[Particle(particle_mass,particle_charge,particle_agg) for p in range((int)(particle_num/particle_agg))]


    def create_mesh(self,delta_x,num_cells):
        self.mesh = [Node(n*self.delta_x) for n in range(self.num_cells)]



    def create_poisson_matrix(self):
        A=np.zeros((n,m))
        

    