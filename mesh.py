import numpy as np

class Mesh:
    
    def __init__(self, cell_length, cell_number):
        self.delta_x=cell_length
        self.charges=np.zeros((cell_number,cell_number))
        self.fields_x=np.zeros((cell_number,cell_number))
        self.fields_y=np.zeros((cell_number,cell_number))

    def zero_charge(self):
        self.charges*=0
    
    def zero_field(self):
        self.fields_x*=0
        self.fields_y*=0
    
    def aggregate_charge(self, bottom_left, particles): # surrounding_* is an array of 4-tuples, in order of top left, top right, bottom left, bottom right
        for i in range(particles.number):
           pass 



            # self.charges[left_nodes[i]]+=particles.charges[i] * left_weights[i]
            # self.charges[right_nodes[i]]+=particles.charges[i] * right_weights[i]

