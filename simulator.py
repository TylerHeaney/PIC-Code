import numpy as np
from particles import ParticleHandler
from mesh import Mesh
import random
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve


epsilon_naught=8.85418e-12

class Simulator:
    """The actual simulator"""

    def shapes(self,distances):
       pass 

    def __init__(self,delta_t,delta_x,num_cells,particle_num,particle_mass,particle_charge,particle_agg):
        self.step_num=0
        self.delta_t=delta_t
        self.delta_x=delta_x
        self.num_cells=num_cells
        self.particles=ParticleHandler(particle_mass,particle_charge,particle_agg)
        self.mesh=Mesh(delta_x,num_cells)
        self.particle_num=particle_num

        self.A,self.FD = self.create_matrices()


    def create_matrices(self):
        e=np.ones(self.num_cells)

        diags = np.array([-1,0,1])
        vals  = np.vstack((e,-2*e,e))
        A = sp.spdiags(vals, diags, self.num_cells, self.num_cells)
        A = sp.lil_matrix(A)
        A[0,self.num_cells-1] = 1
        A[self.num_cells-1,0] = 1
        A /= self.delta_x**2
        A = sp.csr_matrix(A)

        diags = np.array([-1,1])
        vals = np.vstack((-1*e,e))
        FD = sp.spdiags(vals,diags,self.num_cells,self.num_cells)
        FD = sp.lil_matrix(FD)
        FD[0,self.num_cells-1]=-1
        FD[self.num_cells-1,0]=1
        FD /= (2*self.delta_x)
        FD = sp.csr_matrix(FD)

        return A,FD


    def initialize(self, positions, momenta):
        if(len(positions)!=self.particle_num or len(momenta)!=self.particle_num): # use random numbers
            self.particles.initialize([[random.random()*self.delta_x*self.num_cells]*2 for i in range(self.particle_num)],0)
        else:
            self.particles.initialize(positions,momenta)
        



    def step(self):

        left_col=np.floor(self.particles.positions[:,0]/self.delta_x).astype(int)
        bottom_row=np.floor(self.particles.positions[:,1]/self.delta_x).astype(int)
        bottom_left=np.column_stack((left_col,bottom_row))
        # bilinear interpolation next

        # left_nodes=np.floor(self.particles.positions/self.delta_x).astype(int)
        # right_nodes=left_nodes+1
        # left_weights=self.shapes(left_nodes*self.delta_x-self.particles.positions)
        # right_weights=self.shapes(right_nodes*self.delta_x-self.particles.positions)
        # right_nodes=np.mod(right_nodes,self.num_cells)


        self.scatter(bottom_left)
        # self.field_solve()
        self.gather(bottom_left)
        # self.push()
        self.step_num+=1
   
   
    def scatter(self,bottom_left):
        self.mesh.zero_charge()
        self.mesh.aggregate_charge(bottom_left,self.particles)
    

    def field_solve(self):
        self.mesh.zero_field()
        potential_vector=spsolve(self.A,self.charge_vector()-np.sum(self.charge_vector())/self.num_cells,permc_spec="MMD_AT_PLUS_A")
        self.find_E(potential_vector)


    def charge_vector(self):
        b=-1*self.mesh.charges/epsilon_naught
        return b
    

    def find_E(self, potential_vector):

        self.mesh.fields = -1 * self.FD @ potential_vector
        
        


    def gather(self,bottom_left):
        self.particles.zero_field()
        self.particles.aggregate_field(bottom_left, self.mesh, self.delta_x)


    def push(self):
        self.particles.momenta+=self.particles.charges*self.particles.fields*self.delta_t
        self.particles.positions+=self.particles.velocities()*self.delta_t
        self.particles.positions=np.mod(self.particles.positions, self.num_cells*self.delta_x)
        




    