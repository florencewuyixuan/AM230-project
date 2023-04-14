# %matplotlib widget
# import gpumd as md
import time
import json

from diffractio import sp, nm, plt, np, mm, degrees, um
from diffractio.scalar_fields_XY import Scalar_field_XY
from diffractio.scalar_sources_XY import Scalar_source_XY
from diffractio.scalar_masks_XY import Scalar_mask_XY
from diffractio.utils_drawing import draw_several_fields

from random import uniform
from math import pi, sin, cos

import matplotlib.pyplot as plt
import pdb

Nsteps = 1000
if 0:
    for t in range(Nsteps):
        json_reader = md.fast_read_json("initphi=0.4L=150.json") #here we read the json file in C++
        system = md.System(json_reader.particles, json_reader.box)

        # #get the GPU hardware properties
        # cuda_info = system.get_execution_policies()
        # print("** Device parameters **")
        # print(cuda_info.getDeviceProperty())
        # print("** Lunch parameters **")
        # print("GridSize=", cuda_info.getGridSize(), "BlockSize=", cuda_info.getBlockSize())

        dump = md.Dump(system)          # Create a dump object
        print("t=0")
        dump.show(notebook=True)        # Plot the particles with matplotlib

        evolver = md.Evolver(system)    # Create a system evolver object

        #add the forces and torques

        # Create pairwise repulsive interactions with the spring contant k = 10 and range a = 1.0
        evolver.add_force("Harmonic Force", {"k": 10.0, "a": 1.0})

        # Create self-propulsion, self-propulsion strength alpha = 10.0
        evolver.add_force("Self Propulsion", {"alpha": 10.0})

        # Create pairwise polar alignment with alignment strength J = 10.0 and range a = 1.0
        evolver.add_torque("Polar Align", {"k": 10.0, "a": 1.0})

        #add integrators
        # Integrator for updating particle position, friction gamma = 1.0 , "random seed" seed = 10203 and no thermal noise
        evolver.add_integrator("Brownian Positions", {"T": 0.0, "gamma": 1.0, "seed": 10203})

        # Integrator for updating particle orientation, friction gamma = 1.0, "rotation" T = 0.1, D_r = 0.0, "random seed" seed = 10203
        evolver.add_integrator("Brownian Rotation", {"T": 0.1, "gamma": 1.0, "seed": 10203})

        evolver.set_time_step(1e-2) # Set the time step for all the integrators

        evolver.evolve()    # Evolve the system by one time step

        dump.dump_json(filename)
        
        mask = compute_mask(filename) # Dima
        field = diffractio.propagate(mask) # Dima
        velocity = convertTovVelocity(field) # Dima

        update_json(filename, velocity) # Florence

if 1:
    

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

        num_data = 2**8
        length = 60 * um
        x = np.linspace(-length / 2, length / 2, num_data)
        y = np.linspace(-length / 2, length / 2, num_data)
        wavelength = 0.6328 * um

        

        # pdb.set_trace()

        t_total = Scalar_mask_XY(x, y, wavelength)
        for ii in range(0,len(particles_x[0:100])-1):
            t1 = Scalar_mask_XY(x, y, wavelength)
            t1.circle(r0=(particles_x[ii] * um, particles_y[ii] * um), radius=2*um)
            # t1.draw(kind='intensity')

            t_total = t_total.add(t1,'maximum1')
        
        plt.figure()
        t_total.draw(title = 'Particle mask')
        plt.show()

        field = t_total.RS(z = 0.1 * mm,new_field=True)
        plt.figure()
        field.draw()
        plt.show()

    filename = 'test.json'
    compute_mask(filename)

# diffractio tutorial adding masks
if 0:
    num_data = 512
    length = 250 * um
    x = np.linspace(-length / 2, length / 2, num_data)
    y = np.linspace(-length / 2, length / 2, num_data)
    wavelength = 0.6328 * um

    t1 = Scalar_mask_XY(x, y, wavelength)
    t1.square(r0=(-50 * um, 0 * um), size=(50 * um, 50 * um), angle=0 * degrees)

    t2 = Scalar_mask_XY(x, y, wavelength)
    t2.circle(r0=(50 * um, 0 * um), radius=(25 * um, 25 * um), angle=0 * degrees)

    t3 = t2 + t1

    plt.figure()
    draw_several_fields([t1, t2, t3], titles=['1', '2', '1+2'])
    plt.show()