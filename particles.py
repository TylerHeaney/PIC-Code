import numpy

class Particle:
    """Class for each particle in the simulation"""
    
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
    
    def aggregate_field(self,mesh):
        for node in mesh:
            self.field+=node.field * shape(node.position-self.position)
    
    def velocity(self):
        return momentum/mass


