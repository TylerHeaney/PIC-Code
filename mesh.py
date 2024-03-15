import numpy as np

class Mesh:
    
    def __init__(self, cell_length, cell_number):
        self.delta_x=cell_length
        self.charges=np.zeros(cell_number)
        self.fields=np.zeros(cell_number)

    def zero_charge(self):
        self.charges*=0
    
    def zero_field(self):
        self.fields*=0
    
    def aggregate_charge(self, left_nodes, right_nodes, left_weights, right_weights, particles):
        for i in range(particles.number):
            self.charges[left_nodes[i]]+=particles.charges[i] * left_weights[i]
            self.charges[right_nodes[i]]+=particles.charges[i] * right_weights[i]

