from Python.pymd.builder import *
from Python.pymd.md import *               # Import the md module from the pymd package
import numpy as np
from diffractio import sp, nm, plt, np, mm, degrees, um
from diffractio.scalar_fields_XY import Scalar_field_XY
from diffractio.scalar_sources_XY import Scalar_source_XY
from diffractio.scalar_masks_XY import Scalar_mask_XY
from diffractio.utils_drawing import draw_several_fields
from compute_mask import *

from random import uniform
from math import pi, sin, cos

import matplotlib.pyplot as plt
import pdb

# initializing
phi = 0.4
L = 35
a = 1.0
random_init(phi, L, rcut=a, outfile='testing_0.json')  

# number of iteration
N = 1000

# TODO: loop starts here
for i in range(1, N+1):

    s = System(rcut = 3.0, pad = 0.5)   # Create a system object with neighbour list cutoff rcut = 3.0 and padding distance 0.5
    s.read_init('testing_' + str(i-1) + '.json')            # Read in the initial configuration

    e = Evolver(s)                      # Create a system evolver object
    d = Dump(s)                         # Create a dump object

    hf = HarmonicForce(s, 10.0, 2.0)    # Create pairwise repulsive interactions with the spring contant k = 10 and range a = 2.0

    # TODO:: change alpha to intensity
    # alpha = compute_mask('testing.json')
    alpha = compute_mask('testing_' + str(i-1) + '.json')

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

    # evolve once with dt, dt = 0.01 now
    dt = 0.01
    e.evolve(dt)

    # dump the particles as a json file after each iteration
    d.dump_json('testing_' + str(i) + '.json')
    d.dump_vtp('test_{:05d}.vtp'.format(i))
    ## TODO: we may want to dump to vtk for plot after some iterations

