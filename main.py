import time
from random import random
from math import sqrt
import multiprocessing


class Molecule:
    def __init__(self, radius, V, planeSize):
        self.Radius = radius
        self.randomizePosition(planeSize)
        self.randomizeVelocity(V, planeSize[2])

    def randomizePosition(self, planeSize):
        x, y, z = random() * planeSize[0], random() * planeSize[1], random() * planeSize[2]
        self.Position = (x, y, z)

    def randomizeVelocity(self, V, Z):
        x = random() * V * (1 if random() > 0.5 else -1)
        y = random() * sqrt(V ** 2 - x ** 2) * (1 if random() > 0.5 else -1) if Z else sqrt(V ** 2 - x ** 2) * (1 if random() > 0.5 else -1)
        z = sqrt(V ** 2 - x ** 2 - y ** 2) * (1 if random() > 0.5 else -1) if Z else 0

        self.Velocity = (x, y, z)

    def velocityMod(self):
        return sqrt(self.Velocity[0] ** 2 + self.Velocity[1] ** 2 + self.Velocity[2] ** 2 )

    def distance(self, other):
        return sqrt((self.Position[0] - other.Position[0]) ** 2 + (self.Position[1] - other.Position[1]) ** 2 + (self.Position[2] - other.Position[2]) ** 2 )

    def update(self, dt):
        self.Position = (self.Position[0] + self.Velocity[0] * dt, self.Position[1] + self.Velocity[1] * dt, self.Position[2] + self.Velocity[2] * dt)

    def collision(self, other):
        dist = self.distance(other)
        if dist == 0:
            dist = 0.0000000000000000000000001
        if dist <= self.Radius * 2:
            dx, dy, dz = self.Position[0] - other.Position[0], self.Position[1] - other.Position[1], self.Position[2] - other.Position[2]
            normal = (dx/dist, dy/dist, dz/dist)
            overlap = (dist - (self.Radius * 2)) / 2
            
            self.Position, other.Position = tuple([pos - nor*overlap for pos, nor in zip(self.Position, normal)]), tuple([pos + nor*overlap for pos, nor in zip(other.Position, normal)])
            relVel = tuple([vel[0]-vel[1] for vel in zip(self.Velocity, other.Velocity)])
            RVnMOD = sum([vel*norm for vel, norm in zip(relVel, normal)])
            RVn    = tuple(vel*RVnMOD for vel in normal)
                                            
            self.Velocity, other.Velocity = tuple([vel-rel for vel, rel in zip(self.Velocity, RVn)]), tuple([vel+rel for vel, rel in zip(other.Velocity, RVn)])

    def __str__(self):
        return f'Molecule with position of {self.Position}, velocity vector of {self.Velocity} and velocity module of {self.velocityMod()}'


class Box:
    def __init__(self, size):
        self.Size = size

    def collision(self, molecule, sim):
    
        no3D = (self.Size[2]==0) 
        D01, D12, D02 = molecule.Velocity[0] / molecule.Velocity[1], 0 if no3D else molecule.Velocity[1] / molecule.Velocity[2], 0 if no3D else molecule.Velocity[0] / molecule.Velocity[2]
        if molecule.Position[0] <= 0:
            molecule.Position = (
            0, 
            -molecule.Position[0] / D01 * (1 if molecule.Velocity[1] > 0 else -1) + molecule.Position[1], 
            0 if no3D else -molecule.Position[0] / D02 * (1 if molecule.Velocity[2] > 0 else -1) + molecule.Position[2])
            
            molecule.Velocity = (-molecule.Velocity[0], molecule.Velocity[1],  molecule.Velocity[2])
            
        elif molecule.Position[0] >= self.Size[0]:
            molecule.Position = (
            self.Size[0],
            (molecule.Position[0] - self.Size[0]) / D01 * (1 if molecule.Velocity[1] > 0 else -1) + molecule.Position[1], 
            0 if no3D else (molecule.Position[0] - self.Size[0]) / D02 * (1 if molecule.Velocity[2] > 0 else -1) + molecule.Position[2])
                                 
            molecule.Velocity = (-molecule.Velocity[0], molecule.Velocity[1], molecule.Velocity[2])
            
        elif molecule.Position[1] <= 0:
            molecule.Position = (
            -molecule.Position[1] * D01 * (1 if molecule.Velocity[0] > 0 else -1) + molecule.Position[0], 
            0, 
            0 if no3D else -molecule.Position[1] / D12 * (1 if molecule.Velocity[2] > 0 else -1) + molecule.Position[2])
            
            molecule.Velocity = (molecule.Velocity[0], -molecule.Velocity[1], molecule.Velocity[2])
            
        elif molecule.Position[1] >= self.Size[1]:
            molecule.Position = (
            (molecule.Position[1] - self.Size[1]) * D01 * (1 if molecule.Velocity[0] > 0 else -1) + molecule.Position[0],
            self.Size[1], 
            0 if no3D else (molecule.Position[1] - self.Size[1]) / D12 * (1 if molecule.Velocity[2] > 0 else -1) + molecule.Position[2])
            
            molecule.Velocity = (molecule.Velocity[0], -molecule.Velocity[1], molecule.Velocity[2])
            
        if no3D:
            return
            
        if molecule.Position[2] <= 0:
            molecule.Position = (
            -molecule.Position[2] * D02 * (1 if molecule.Velocity[0] > 0 else -1) + molecule.Position[0], 
            -molecule.Position[2] * D12 * (1 if molecule.Velocity[1] > 0 else -1) + molecule.Position[1], 
            0)            
            
            molecule.Velocity = (molecule.Velocity[0], molecule.Velocity[1], -molecule.Velocity[2])
            
        elif molecule.Position[2] >= self.Size[2]:
            molecule.Position = (
            (molecule.Position[2] - self.Size[2]) * D02 * (1 if molecule.Velocity[0] > 0 else -1) + molecule.Position[0],
            (molecule.Position[2] - self.Size[2]) * D12 * (1 if molecule.Velocity[1] > 0 else -1) + molecule.Position[1], 
            self.Size[2])
            
            molecule.Velocity = (molecule.Velocity[0], molecule.Velocity[1], -molecule.Velocity[2])
            
        else:
            return
        
            
            
        if CAN_LEAVE and molecule.velocityMod() > VELOCITY_LIMIT:
            sim.molecules.remove(molecule)
            del molecule


class Simulation:

    def __init__(self, num, V, radius, planeSize):
        self.dt = radius/V*0.25
        self.N = num
        self.molecules = self.generateFirstState(num, V, radius, planeSize)
        self.box = Box(planeSize)
        self.time_passed = 0
        processes = []
        p_run = multiprocessing.Process(target=self.run)
        p_run.start()
        processes.append(p_run)
        if GRAPH:
            p_graph = multiprocessing.Process(target=graph, args=(VELOCITY_LIMIT,))
            p_graph.start()
            processes.append(p_graph)
        if VISUALISE:
            p_vis = multiprocessing.Process(target=visualize, args=(BOX, RADIUS))
            p_vis.start()
            processes.append(p_vis)
        for p in processes:
            p.join()

    def calc_energy(self):
        E = 0
        for mol in self.molecules:
            E += mol.velocityMod() ** 2
        return E

    def run(self):
        while 1:
            self.update(self.dt)
            print(self.calc_energy())

    def update(self, dt):
        for molecule in self.molecules:
            molecule.update(dt)
        for molecule in self.molecules:
            self.box.collision(molecule, self)
        for i in range(len(self.molecules)):
            for j in range(i + 1, len(self.molecules)):
                self.molecules[i].collision(self.molecules[j])
        self.time_passed += dt
        data = open('data.txt', 'r+')
        data.close()

        with open(r'data.txt', 'w') as data:
            data.write('\n'.join([str(self.time_passed*1000)]+
                list(map(str, [f'{mol.velocityMod()} {mol.Position[0]} {mol.Position[1]} ' for mol in self.molecules]))))

    def generateFirstState(self, num, V, radius, planeSize):
        return list([Molecule(radius, V, planeSize) for i in range(num)])


NUM             = 1000                  # Number of particles
VELOCITY        = 2000                  # Starting velocity module
RADIUS          = 2                     # Radius of particles
BOX             = (150, 150,  150)      # "Crystal" size
VELOCITY_LIMIT  = 3000                  # The speed at which a molecule can leave a crystallite
CAN_LEAVE       = False                 # True if molecule can leave crystallite, False if not
VISUALISE       = False                 # True if you want to visualize the simulation, False if not. DO NOT SET TO TRUE ON LARGE NUMBERS. Only shows 2D
GRAPH           = True                  # True if you want to see animated histogram of particle velocities, False if not

if __name__ == '__main__':
    from render_sym import visualize
    from graph_data import graph

    sim = Simulation(NUM, VELOCITY, RADIUS, BOX)
