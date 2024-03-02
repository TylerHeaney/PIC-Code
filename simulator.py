import numpy as np
from particles import Particle
from mesh import Node
import random


epsilon_naught=8.85418e-12

class Simulator:
    """The actual simulator"""

    def __init__(self,delta_t,delta_x,num_cells,particle_num,particle_mass,particle_charge,particle_agg):
        self.step_num=0
        self.delta_t=delta_t
        self.delta_x=delta_x
        self.num_cells=num_cells
        self.create_particles(particle_num,particle_agg,particle_mass,particle_charge)
        self.create_mesh(delta_x,num_cells)


    def create_particles(self,particle_num,particle_agg,particle_mass,particle_charge):
        self.particles=[Particle(particle_mass[p],particle_charge[p],particle_agg[p]) for p in range(particle_num)]


    def create_mesh(self,delta_x,num_cells):
        self.mesh = [Node(n*self.delta_x) for n in range(self.num_cells)]
        

    def initialize(self, lhw, rhw, positions_and_momentums):
        if(len(positions_and_momentums)!=len(self.particles)): # use random numbers
            for p in self.particles:
                p.initialize(random.random()*self.delta_x*self.num_cells/2+.25*self.delta_x*self.num_cells,0)
        else:
            for i in range(len(self.particles)):
                self.particles[i].initialize(*positions_and_momentums[i])
        
        self.lhw=lhw
        self.mesh[0].charge=lhw
        self.rhw=rhw
        self.mesh[-1].charge=rhw

    def step(self):
        self.scatter()
        self.field_solve()
        self.gather()
        self.push()
        self.step_num+=1
   
   
    def scatter(self):
        self.zero_charge()
        for node in self.mesh:
            node.aggregate_charge(self.particles,self.delta_x)
    

    def zero_charge(self):
        for node in self.mesh:
            node.charge=0


    def field_solve(self):
        self.mesh_zero_field()
        A=self.create_poisson_matrix()
        potential_vector=np.linalg.solve(A,self.charge_vector())
        self.find_E(potential_vector)


    def mesh_zero_field(self):
        for node in self.mesh:
            node.field=0


    def create_poisson_matrix(self):
        A=np.zeros((len(self.mesh),len(self.mesh)))
        A[0,0]=1
        for i in range(1,len(self.mesh)-1):
            A[i,i-1]=1/self.delta_x**2
            A[i,i]=-2/self.delta_x**2
            A[i,i+1]=1/self.delta_x**2
        A[-1,-1]=1;
        return A


    def charge_vector(self):
        b=np.zeros(len(self.mesh))
        for i in range(len(self.mesh)):
            b[i]=-self.mesh[i].charge/epsilon_naught
        return b
    

    def find_E(self, potential_vector):
        self.mesh[0].field=self.lhw
        for i in range(1,len(self.mesh)-1):
            self.mesh[i].field=-1*(potential_vector[i+1]-potential_vector[i-1])/2/self.delta_x
        self.mesh[-1].field=self.rhw


    def gather(self):
        self.zero_field()
        for particle in self.particles:
            particle.aggregate_field(self.mesh,self.delta_x)
    

    def zero_field(self):
        for p in self.particles:
            p.field=0


    def push(self):
        for particle in self.particles:
            if particle.position < 0 or particle.position > self.num_cells*self.delta_x:
                particle.momentum=-1*particle.momentum
            particle.momentum=particle.momentum+particle.charge*particle.field*self.delta_t
            particle.position=particle.position+particle.velocity()*self.delta_t

        




    