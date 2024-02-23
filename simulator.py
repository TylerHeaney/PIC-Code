import numpy as np
import particles
import mesh

class Simulator:
    """The actual simulator"""

    def __init__(self,delta_t,delta_x,cell_num,rhw,lhw,particle_num,particle_mass,particle_charge,particle_agg,positions,momentums):
        self.step_num=0
        self.delta_t=delta_t
        self.delta_x=delta_x
        self.cell_num=cell_num
        self.particles = self.create_particles(particle_num,particle_agg,particle_mass,particle_charge)
        self.initialize_particles(positions,momentums)
        self.mesh = self.create_mesh(delta_x,cell_num)
        self.initialize_mesh(rhw,lhw)
    
    def step(self):
        self.scatter()
        self.field_solve()
        self.gather()
        self.push()


        self.step_num+=1
   
   
    def scatter(self):
        mesh.zero_charge()
        for node in mesh:
            node.aggregate_charge(particles)
    
    
    def field_solve(self):
        mesh.zero_field()
        A=self.create_poisson_matrix
        potential_vector=np.linalg.solve(a,mesh.charge_vector)
        mesh.find_E(potential_vector)
    
    def gather(self):
        particles.zero_field()
        for particle in particles:
            particles.aggregate_field(mesh)
    

    def push(self):
        for particle in particles:
            particle.position=particle.position+particle.velocity()*self.time_step()
            particle.momentum=particle.momentum+particle.charge*particle.field*self.time_step()


    def time_step(self):
        self.step_num*delta_t