import numpy

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


