import numpy as np

class ParticleHandler:

    def __init__(self,masses,charges,agg):
        if len(masses)!=len(charges):
            raise ValueError("masses and charges must be equal")
        self.masses=np.array(masses)*agg
        self.charges=np.array(charges)*agg
        self.agg=agg
        self.positions=np.zeros((len(masses),2))
        self.momenta=np.zeros((len(masses),2))
        self.fields_x=np.zeros((len(masses),2))
        self.fields_y=np.zeros((len(masses),2))
        self.number=len(masses)

    def initialize(self, positions, momenta):
        if(len(positions)!=len(momenta) or len(positions)!=len(self.masses)):
            raise ValueError("positions and momenta must be the same length as number of particles")
        
        self.positions=np.array(positions)
        self.momenta=np.array(momenta)

    def zero_field(self):
        self.fields_x*=0
        self.fields_y*=0

    def aggregate_field(self, lower_left,rel_coord,mesh,delta_x):
        for i in range(len(weights)):
            

            # self.fields[i]+=mesh.fields[left_nodes[i]] * left_weights[i]
            # self.fields[i]+=mesh.fields[right_nodes[i]] * right_weights[i]

    def size(self):
        return len(self.masses)

    def velocities(self):
        return self.momenta/self.masses