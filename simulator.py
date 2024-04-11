import numpy as np
from particles import ParticleHandler
from mesh import Mesh
import random
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve


epsilon_naught=8.85418e-12/100

class Simulator:
    """The actual simulator"""

    def shapes(self,distances):
       return 1-distances

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
        vals  = np.vstack((-1*e,4*e,-1*e))
        D = sp.spdiags(vals, diags, self.num_cells, self.num_cells)
        D = sp.lil_matrix(D)
        D[0,self.num_cells-1] = -1
        D[self.num_cells-1,0] = -1
        # D /= self.delta_x**2
        D = sp.csr_matrix(D).toarray()

        I=-1*np.identity(self.num_cells)
        A=np.zeros([self.num_cells*self.num_cells]*2)
        # A=np.bmat([[I]*self.num_cells]*self.num_cells)
        for i in range(1,self.num_cells-1):
            A[i*self.num_cells:i*self.num_cells+self.num_cells,i*self.num_cells:i*self.num_cells+self.num_cells]=D
            A[(i-1)*self.num_cells:(i-1)*self.num_cells+self.num_cells,i*self.num_cells:i*self.num_cells+self.num_cells]=I
            A[(i+1)*self.num_cells:(i+1)*self.num_cells+self.num_cells,i*self.num_cells:i*self.num_cells+self.num_cells]=I
        # A[-1*self.num_cells:,0:self.num_cells]=I
        # A[0:self.num_cells,-1*self.num_cells:]=I
        A[0*self.num_cells:0*self.num_cells+self.num_cells,0*self.num_cells:0*self.num_cells+self.num_cells]=D
        A[-1*self.num_cells:,-1*self.num_cells:]=D

        A=sp.csr_matrix(A)

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
        dx=self.particles.positions[:,0]/self.delta_x-left_col
        dy=self.particles.positions[:,1]/self.delta_x-bottom_row
        blw=self.shapes(dx)*self.shapes(dy)
        brw=self.shapes(1-dx)*self.shapes(dy)
        tlw=self.shapes(dx)*self.shapes(1-dy)
        trw=self.shapes(1-dx)*self.shapes(1-dy)
        # bilinear interpolation next

        # left_nodes=np.floor(self.particles.positions/self.delta_x).astype(int)
        # right_nodes=left_nodes+1
        # left_weights=self.shapes(left_nodes*self.delta_x-self.particles.positions)
        # right_weights=self.shapes(right_nodes*self.delta_x-self.particles.positions)
        # right_nodes=np.mod(right_nodes,self.num_cells)


        self.scatter(bottom_left,blw,brw,tlw,trw)
        self.field_solve()
        self.gather(bottom_left,blw,brw,tlw,trw)
        self.push()
        self.step_num+=1
   
   
    def scatter(self,bottom_left,blw,brw,tlw,trw):
        self.mesh.zero_charge()
        self.mesh.aggregate_charge(bottom_left,blw,brw,tlw,trw,self.particles)
    

    def field_solve(self):
        self.mesh.zero_field()
        potential_matrix = np.reshape(spsolve(self.A,self.charge_vector(),permc_spec="MMD_AT_PLUS_A"),(self.num_cells,self.num_cells))
        self.find_E(potential_matrix)



    def charge_vector(self):
        b=-1*self.mesh.charges/epsilon_naught
        return b.flatten()[:,None]
    

    def find_E(self, potential_matrix):
        self.mesh.fields[:,:,0] = (-1 * self.FD @ potential_matrix.T).T
        self.mesh.fields[:,:,1] = -1 * self.FD @ potential_matrix
        
        


    def gather(self,bottom_left,blw,brw,tlw,trw):
        self.particles.zero_field()
        self.particles.aggregate_field(bottom_left,blw,brw,tlw,trw, self.mesh, self.delta_x)


    def push(self):
        self.particles.momenta+=self.particles.charges[:,None]*self.particles.fields*self.delta_t
        self.particles.positions+=self.particles.velocities()*self.delta_t
        self.particles.positions=np.mod(self.particles.positions, self.num_cells*self.delta_x)
        




    