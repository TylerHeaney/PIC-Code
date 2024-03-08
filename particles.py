import numpy as np

class Particle:
    """Class for each particle in the simulation"""
    
    def shape(distance):
        distance=abs(distance)
        if(distance>1):
            return 0
        return distance

    def __init__(self,real_mass,real_charge,number):
        self.mass=real_mass*number
        self.charge=real_charge*number
        self.number=number
        self.position=float('nan')
        self.momentum=float('nan')
        self.field=float('nan')
        
    def initialize(self,position, momentum):
        self.position=position
        self.momentum=momentum
    
    def aggregate_field(self,mesh, delta_x):
        for node in mesh:
            self.field+=node.field * Particle.shape((node.position-self.position)/delta_x)
    
    def velocity(self):
        return self.momentum/self.mass


class ParticleHandler:

    def shape(distance):
        distance=abs(distance)
        if(distance>1):
            return 0
        return distance

    def __init__(self,masses,charges,agg):
        if len(masses)!=len(charges):
            raise ValueError("masses and charges must be equal")
        self.masses=np.array(masses)*agg
        self.charges=np.array(charges)*agg
        self.agg=agg
        self.positions=np.zeros(len(masses))
        self.momenta=np.zeros(len(masses))
        self.fields=np.zeros(len(masses))
        self.number=len(masses)

    def initialize(self, positions, momenta):
        if(len(positions)!=len(momenta) or len(positions)!=len(self.masses)):
            raise ValueError("positions and momenta must be the same length as number of particles")
        
        self.positions=np.array(positions)
        self.momenta=np.array(momenta)

    def zero_field(self):
        self.fields*=0

    def aggregate_field(self,left_nodes,right_nodes,left_weights,right_weights,mesh,delta_x):
        for i in range(len(left_nodes)):
            self.fields[i]+=mesh.fields[left_nodes[i]] * left_weights[i]
            self.fields[i]+=mesh.fields[right_nodes[i]] * right_weights[i]

    def size(self):
        return len(self.masses)

    def velocities(self):
        return self.momenta/self.masses