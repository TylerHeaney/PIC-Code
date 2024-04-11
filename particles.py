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
        self.fields=np.zeros((len(masses),2))
        self.number=len(masses)

    def initialize(self, positions, momenta):
        if(len(positions)!=len(momenta) or len(positions)!=len(self.masses)):
            raise ValueError("positions and momenta must be the same length as number of particles")
        
        self.positions=np.array(positions)
        self.momenta=np.array(momenta)

    def zero_field(self):
        self.fields*=0

    def aggregate_field(self, bottom_left,blw,brw,tlw,trw,mesh,delta_x):
        bottom_right=(bottom_left+[1,0])%mesh.cell_number
        top_left=(bottom_left+[0,1])%mesh.cell_number
        top_right=(bottom_left+[1,1])%mesh.cell_number
        for i in range(len(bottom_left)):
            self.fields[i]+=mesh.fields[bottom_left[i][0],bottom_left[i][1]] * blw[i]
            self.fields[i]+=mesh.fields[(bottom_left[i][0]+1)%mesh.cell_number,bottom_left[i][1]] * brw[i]
            self.fields[i]+=mesh.fields[bottom_left[i][0],(bottom_left[i][1]+1)%mesh.cell_number] * tlw[i]
            self.fields[i]+=mesh.fields[(bottom_left[i][0]+1)%mesh.cell_number,(bottom_left[i][1]+1)%mesh.cell_number] * trw[i]
            

            # self.fields[i]+=mesh.fields[left_nodes[i]] * left_weights[i]
            # self.fields[i]+=mesh.fields[right_nodes[i]] * right_weights[i]

    def size(self):
        return len(self.masses)

    def velocities(self):
        return self.momenta/[[mass,mass] for mass in self.masses]