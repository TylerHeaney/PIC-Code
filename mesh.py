import numpy as np

class Mesh:
    
    def __init__(self, cell_length, cell_number):
        self.delta_x=cell_length
        self.charges=np.zeros((cell_number,cell_number))
        self.fields=np.zeros((cell_number,cell_number,2))
        self.cell_number=cell_number

    def zero_charge(self):
        self.charges*=0
    
    def zero_field(self):
        self.fields*=0
    
    def aggregate_charge(self, bottom_left,blw,brw,tlw,trw, particles): 
        bottom_right=(bottom_left+[1,0])%self.cell_number
        top_left=(bottom_left+[0,1])%self.cell_number
        top_right=(bottom_left+[1,1])%self.cell_number
        for i in range(particles.number):
            self.charges[bottom_left[i][0],bottom_left[i][1]]+=particles.charges[i]*blw[i]
            self.charges[(bottom_left[i][0]+1)%self.cell_number,bottom_left[i][1]]+=particles.charges[i]*brw[i]
            self.charges[bottom_left[i][0],(bottom_left[i][1]+1)%self.cell_number]+=particles.charges[i]*tlw[i]
            self.charges[(bottom_left[i][0]+1)%self.cell_number,(bottom_left[i][1]+1)%self.cell_number]+=particles.charges[i]*trw[i]




            # self.charges[left_nodes[i]]+=particles.charges[i] * left_weights[i]
            # self.charges[right_nodes[i]]+=particles.charges[i] * right_weights[i]

