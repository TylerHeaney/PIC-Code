import numpy as np
from particles import ParticleHandler
from mesh import Mesh
import random


epsilon_naught=8.85418e-12

class Simulator:
    """The actual simulator"""

    def shapes(self,distances):
        return np.abs(distances)/self.delta_x

    def __init__(self,delta_t,delta_x,num_cells,particle_num,particle_mass,particle_charge,particle_agg):
        self.step_num=0
        self.delta_t=delta_t
        self.delta_x=delta_x
        self.num_cells=num_cells
        self.particles=ParticleHandler(particle_mass,particle_charge,particle_agg)
        self.mesh=Mesh(delta_x,num_cells)
        self.particle_num=particle_num

        self.A = self.create_poisson_matrix()


    def create_poisson_matrix(self):
        A=np.zeros((self.num_cells,self.num_cells))
        A[0,0]=-2/self.delta_x**2
        A[0,1]=1/self.delta_x**2
        A[0,-1]=1/self.delta_x**2
        for i in range(1,self.num_cells-1):
            A[i,i-1]=1/self.delta_x**2
            A[i,i]=-2/self.delta_x**2
            A[i,i+1]=1/self.delta_x**2
        A[-1,-1]=-2/self.delta_x**2
        A[-1,0]=1/self.delta_x**2
        A[-1,-2]=1/self.delta_x**2
        return A


    def initialize(self, lhw, rhw, positions, momenta):
        if(len(positions)!=self.particle_num or len(momenta)!=self.particle_num): # use random numbers
            self.particles.initialize([random.random()*self.delta_x*self.num_cells for i in range(self.particle_num)],0)
        else:
            self.particles.initialize(positions,momenta)
        
        self.lhw=lhw
        self.mesh.charges[0]=lhw
        self.rhw=rhw
        self.mesh.charges[-1]=rhw


    def step(self):
        left_nodes=np.floor(self.particles.positions/self.delta_x).astype(int)
        right_nodes=left_nodes+1
        left_weights=self.shapes(left_nodes*self.delta_x-self.particles.positions)
        right_weights=self.shapes(right_nodes*self.delta_x-self.particles.positions)
        right_nodes=np.mod(right_nodes,self.num_cells)
        self.scatter(left_nodes,right_nodes,left_weights,right_weights)
        self.field_solve()
        self.gather(left_nodes,right_nodes,left_weights,right_weights)
        self.push()
        self.step_num+=1
   
   
    def scatter(self,left_nodes,right_nodes,left_weights,right_weights):
        self.mesh.zero_charge()
        self.mesh.aggregate_charge(left_nodes,right_nodes,left_weights,right_weights,self.particles)
    

    def field_solve(self):
        self.mesh.zero_field()
        potential_vector=np.linalg.solve(self.A,self.charge_vector())
        self.find_E(potential_vector)


    def charge_vector(self):
        b=-1*self.mesh.charges/epsilon_naught
        return b
    

    def find_E(self, potential_vector):
        self.mesh.fields[0]=-1*(potential_vector[1]-potential_vector[-1])/2/self.delta_x
        for i in range(1,self.num_cells-1):
            self.mesh.fields[i]=-1*(potential_vector[i+1]-potential_vector[i-1])/2/self.delta_x
        self.mesh.fields[-1]=-1*(potential_vector[0]-potential_vector[-2])/2/self.delta_x


    def gather(self,left_nodes,right_nodes,left_weights,right_weights):
        self.particles.zero_field()
        self.particles.aggregate_field(left_nodes,right_nodes,left_weights,right_weights, self.mesh, self.delta_x)


    def push(self):
        self.particles.momenta+=self.particles.charges*self.particles.fields*self.delta_t
        self.particles.positions+=self.particles.velocities()*self.delta_t
        self.particles.positions=np.mod(self.particles.positions, self.num_cells*self.delta_x)
        




    