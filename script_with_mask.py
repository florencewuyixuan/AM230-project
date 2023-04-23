from Python.pymd.builder import *
from Python.pymd.md import *               # Import the md module from the pymd package
import numpy as np
from diffractio import sp, nm, plt, np, mm, degrees, um
from diffractio.scalar_fields_XY import Scalar_field_XY
from diffractio.scalar_sources_XY import Scalar_source_XY
from diffractio.scalar_masks_XY import Scalar_mask_XY
from diffractio.utils_drawing import draw_several_fields

from random import uniform
from math import pi, sin, cos

import matplotlib.pyplot as plt
import pdb

# initializing
phi = 0.4
L = 10
a = 1.0
random_init(phi, L, rcut=a, outfile='testing.json')  


# TODO: loop starts here

s = System(rcut = 3.0, pad = 0.5)   # Create a system object with neighbour list cutoff rcut = 3.0 and padding distance 0.5
s.read_init('init.json')            # Read in the initial configuration

e = Evolver(s)                      # Create a system evolver object
d = Dump(s)                         # Create a dump object

hf = HarmonicForce(s, 10.0, 2.0)    # Create pairwise repulsive interactions with the spring contant k = 10 and range a = 2.0

# TODO:: change alpha to intensity
# alpha = compute_mask('testing.json')
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

# evolve once with dt, dt = 1 now
e.evolve(1)

# TODO: dump the particles as a json file after each iteration
# d.dump_json(outfile)


# compute mask
def compute_mask(filename):
        with open(filename) as f:
            data = json.load(f)
            if 'box' in data['system']:
                Lx = data["system"]["box"]["Lx"]
                Ly = data["system"]["box"]["Ly"]
                #   self.box = Box(Lx, Ly)
            else:
                raise Exception('Input JSON file has to include system box section.')
            if 'particles' in data['system']:
                particles_x = []
                particles_y = []

            for p in data['system']['particles']:
                idx = p['id']
                x, y = p['r']
                theta = uniform(-pi,pi)
                nx, ny = cos(theta), sin(theta)
                vx, vy = 0.0, 0.0
                fx, fy = 0.0, 0.0
                if 'n' in p:  nx, ny = p['n']
                if 'v' in p:  vx, vy = p['v']
                if 'f' in p:  fx, fy = p['f']
                particles_x.append(x)
                particles_y.append(y)

        num_data = 2**9
        length = 60 * um
        x = np.linspace(-length / 2, length / 2, num_data)
        y = np.linspace(-length / 2, length / 2, num_data)
        wavelength = 0.6328 * um
        # wavelength = 1.55 * um

        

        # pdb.set_trace()
        # N_particles = 200

        # generate scalar mask based on particle positions
        t_total = Scalar_mask_XY(x, y, wavelength)
        #  container to save individual particle masks
        t_list = []
        for ii in range(0,len(particles_x)-1):
            t1 = Scalar_mask_XY(x, y, wavelength)
            t1.circle(r0=(particles_x[ii] * um, particles_y[ii] * um), radius=1*um)
            # t1.draw(kind='intensity')
            t_list.append(t1)

            t_total = t_total.add(t1,'maximum1')
        
        t_total.u = np.abs(1-t_total.u)
        # plt.figure()
        # t_total.draw(title = 'Particle mask')
        # plt.scatter(particles_x[0:N_particles-1],particles_y[0:N_particles-1])
        # plt.show()

        # compute the field distribution after propagating a set distance using Rayleigh-Sommerfeld model
        u0 = Scalar_source_XY(x=x,y=y,wavelength=wavelength)
        u0.gauss_beam([0,0],[100*um,100*um],[0,0])
        # u0.plane_wave(A=1)
        field = (t_total * u0).RS(z = 100 * um,new_field=True)
        
        # field.draw(title = 'Field intensity')
        # cbar = plt.colorbar()
        # plt.scatter(particles_x[0:N_particles-1],particles_y[0:N_particles-1])
        # plt.show()

        
            

        # compute force from intensity gradients
        intensity = np.square(np.abs(field.u))
        # force_x = np.gradient(intensity)[0]
        # force_y = np.gradient(intensity)[1]
        # pdb.set_trace()

        v_update = []
        for p in t_list:
            # integrate intensity within a circle for each particle
            v = np.multiply(intensity,p.u).sum() 
            v_update.append(v)
        
        

        return v_update