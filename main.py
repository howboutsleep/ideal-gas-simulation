import time
from random import random
from math import sqrt
import multiprocessing


class Molecule:
    def __init__(self, radius, V, planeSize):
        self.Radius = radius
        self.randomizePosition(planeSize)
        self.randomizeVelocity(V)

    def randomizePosition(self, planeSize):
        x, y = random() * planeSize[0], random() * planeSize[1]
        self.Position = (x, y)

    def randomizeVelocity(self, V):
        x = random() * V * (1 if random() > 0.5 else -1)
        y = sqrt(V ** 2 - x ** 2) * (1 if random() > 0.5 else -1)
        self.Velocity = (x, y)

    def velocityMod(self):
        return sqrt(self.Velocity[0] ** 2 + self.Velocity[1] ** 2)

    def distance(self, other):
        return sqrt((self.Position[0] - other.Position[0]) ** 2 + (self.Position[1] - other.Position[1]) ** 2)

    def update(self, dt):
        self.Position = (self.Position[0] + self.Velocity[0] * dt, self.Position[1] + self.Velocity[1] * dt)

    def collision(self, other):
        dist = self.distance(other)
        if dist <= self.Radius * 2:
            dx, dy = self.Position[0] - other.Position[0], self.Position[1] - other.Position[1]
            dvx, dvy = self.Velocity[0] - other.Velocity[0], self.Velocity[1] - other.Velocity[1]
            sin, cos = dx / dist, dy / dist
            overlap = (dist - (self.Radius * 2)) / 2
            self.Position, other.Position = (self.Position[0] - sin * overlap, self.Position[1] - cos * overlap), (
            other.Position[0] + sin * overlap, other.Position[1] + cos * overlap)
            self.Velocity, other.Velocity = (self.Velocity[0] - (sin ** 2 * dvx + sin * cos * dvy),
                                             self.Velocity[1] - (sin * cos * dvx + cos ** 2 * dvy)), (
                                            other.Velocity[0] + (sin ** 2 * dvx + sin * cos * dvy),
                                            other.Velocity[1] + (sin * cos * dvx + cos ** 2 * dvy))

    def __str__(self):
        return f'Molecule with position of {self.Position}, velocity vector of {self.Velocity} and velocity module of {self.velocityMod()}'


class Box:
    def __init__(self, size):
        self.Size = size

    def collision(self, molecule, sim):
        tg = molecule.Velocity[0] / molecule.Velocity[1]
        if molecule.Position[0] <= 0:
            molecule.Position = (
            0, -molecule.Position[0] / tg * (1 if molecule.Velocity[1] > 0 else -1) + molecule.Position[1])
            molecule.Velocity = (-molecule.Velocity[0], molecule.Velocity[1])
        elif molecule.Position[0] >= self.Size[0]:
            molecule.Position = (self.Size[0],
                                 (molecule.Position[0] - self.Size[0]) / tg * (1 if molecule.Velocity[1] > 0 else -1) +
                                 molecule.Position[1])
            molecule.Velocity = (-molecule.Velocity[0], molecule.Velocity[1])
        elif molecule.Position[1] <= 0:
            molecule.Position = (
            -molecule.Position[1] * tg * (1 if molecule.Velocity[0] > 0 else -1) + molecule.Position[0], 0)
            molecule.Velocity = (molecule.Velocity[0], -molecule.Velocity[1])
        elif molecule.Position[1] >= self.Size[1]:
            molecule.Position = (
            (molecule.Position[1] - self.Size[1]) * tg * (1 if molecule.Velocity[0] > 0 else -1) + molecule.Position[0],
            self.Size[1])
            molecule.Velocity = (molecule.Velocity[0], -molecule.Velocity[1])
        else:
            return
        if CAN_LEAVE and molecule.velocityMod() > VELOCITY_LIMIT:
            sim.molecules.remove(molecule)
            del molecule


class Simulation:

    def __init__(self, num, V, radius, planeSize):
        self.N = num
        self.molecules = self.generateFirstState(num, V, radius, planeSize)
        self.box = Box(planeSize)
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
        tl = time.perf_counter()
        while 1:
            tr = time.perf_counter()
            self.update(tr - tl)
            tl = tr

    def update(self, dt):
        for molecule in self.molecules:
            molecule.update(dt)
        for molecule in self.molecules:
            self.box.collision(molecule, self)
        for i in range(len(self.molecules)):
            for j in range(i + 1, len(self.molecules)):
                self.molecules[i].collision(self.molecules[j])
        data = open('data.txt', 'r+')
        data.close()

        with open(r'data.txt', 'w') as data:
            data.write('\n'.join(
                map(str, [f'{mol.Position[0]} {mol.Position[1]} {mol.velocityMod()}' for mol in self.molecules])))

    def generateFirstState(self, num, V, radius, planeSize):
        return list([Molecule(radius, V, planeSize) for i in range(num)])


NUM             = 1000         # Number of particles
VELOCITY        = 1000         # Starting velocity module
RADIUS          = 5            # Radius of particles
BOX             = (1000, 1000) # "Crystal" size
VELOCITY_LIMIT  = 2000         # The speed at which a molecule can leave a crystallite
CAN_LEAVE       = False        # True if molecule can leave crystallite, False if not
VISUALISE       = False        # True if you want to visualize the simulation, False if not. DO NOT SET TO TRUE ON LARGE NUMBERS
GRAPH           = True         # True if you want to see animated histogram of particle velocities, False if not

if __name__ == '__main__':
    from render_sym import visualize
    from graph_data import graph

    sim = Simulation(NUM, VELOCITY, RADIUS, BOX)