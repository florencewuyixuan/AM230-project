# import sys
# sys.path.insert(0, '/home/florence/Desktop/AM230/ABPTutorial')
from Python.pymd.builder import *
from Python.pymd.md import *               # Import the md module from the pymd package
import numpy as np


phi = 0.4
L = 10
a = 1.0
random_init(phi, L, rcut=a, outfile='testing.json')  



s = System(rcut = 3.0, pad = 0.5)   # Create a system object with neighbour list cutoff rcut = 3.0 and padding distance 0.5
s.read_init('init.json')            # Read in the initial configuration

e = Evolver(s)                      # Create a system evolver object
d = Dump(s)                         # Create a dump object

hf = HarmonicForce(s, 10.0, 2.0)    # Create pairwise repulsive interactions with the spring contant k = 10 and range a = 2.0
alpha = np.arange(1,41)
sp = SelfPropulsion(s, alpha)         # Create self-propulsion, self-propulsion strength alpha = 1.0
pa = PolarAlign(s, 1.0, 2.0)        # Create pairwise polar alignment with alignment strength J = 1.0 and range a = 2.0

pos_integ = BrownianIntegrator(s, T = 0.0, gamma = 1.0)       # Integrator for updating particle position, friction gamma = 1.0 and no thermal noise
rot_integ = BrownianRotIntegrator(s, T = 0.1, gamma = 1.0)    # Integrator for updating particle oriantation, friction gamma = 1.0, "rotation" T = 0.1, D_r = 0.0 

# Register all forces, torques and integrators with the evolver object
e.add_force(hf)                    
e.add_force(sp)
e.add_torque(pa)
e.add_integrator(pos_integ)
e.add_integrator(rot_integ)

e.evolve(1)

