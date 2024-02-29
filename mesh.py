import numpy

class Node:
    """Discrete points on the Computational Mesh"""

    def shape(distance):
        distance=abs(distance)
        if(distance>1):
            return 0
        return distance

    def __init__(self, position):
        self.position=position
        self.charge=0
        self.field=0

    def aggregate_charge(self, particles):
        for p in particles:
            self.charge+=p.charge * Node.shape(p.position-self.position)

