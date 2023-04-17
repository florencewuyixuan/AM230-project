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

        num_data = 2**9
        length = 60 * um
        x = np.linspace(-length / 2, length / 2, num_data)
        y = np.linspace(-length / 2, length / 2, num_data)
        wavelength = 0.6328 * um
        # wavelength = 1.55 * um

        

        # pdb.set_trace()
        N_particles = 200

        # generate scalar mask based on particle positions
        t_total = Scalar_mask_XY(x, y, wavelength)
        #  container to save individual particle masks
        t_list = []
        for ii in range(0,len(particles_x[0:N_particles])-1):
            t1 = Scalar_mask_XY(x, y, wavelength)
            t1.circle(r0=(particles_x[ii] * um, particles_y[ii] * um), radius=1*um)
            # t1.draw(kind='intensity')
            t_list.append(t1)

            t_total = t_total.add(t1,'maximum1')
        
        t_total.u = np.abs(1-t_total.u)
        # plt.figure()
        t_total.draw(title = 'Particle mask')
        plt.scatter(particles_x[0:N_particles-1],particles_y[0:N_particles-1])
        plt.show()

        # compute the field distribution after propagating a set distance using Rayleigh-Sommerfeld model
        u0 = Scalar_source_XY(x=x,y=y,wavelength=wavelength)
        u0.gauss_beam([0,0],[100*um,100*um],[0,0])
        # u0.plane_wave(A=1)
        field = (t_total * u0).RS(z = 100 * um,new_field=True)
        
        field.draw(title = 'Field intensity')
        cbar = plt.colorbar()
        plt.scatter(particles_x[0:N_particles-1],particles_y[0:N_particles-1])
        plt.show()
        # pdb.set_trace()

        
            

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
        
        # pdb.set_trace()

        # X, Y = np.mgrid[0:len(x),0:len(y)]/len(x)*length - length/2

        # fig, (ax0, ax1, ax2) = plt.subplots(1, 3)
        # ax0.pcolor(X, Y, intensity,cmap=plt.cm.gist_heat)
        # ax1.pcolor(X, Y, force_x,cmap=plt.cm.bone)
        # ax2.pcolor(X, Y, force_y,cmap=plt.cm.bone)

        # plt.show()

        # plt.figure()
        # plt.imshow(x,y,intensity)
        # ax = plt.gca()
        # ax.invert_yaxis()
        # plt.show()

        # plt.figure()
        # plt.subplot(1,2,1)
        # plt.imshow(x,y,force_x)
        # plt.subplot(1,2,2)
        # plt.imshow(x,y,force_y)
        # plt.show()

        # return array of intensities ordered by particle number

        return v_update

    filename = 'test.json'
    v_update = compute_mask(filename)
    pdb.set_trace()
    print()

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

    # plt.figure()
    draw_several_fields([t1, t2, t3], titles=['1', '2', '1+2'])
    plt.show()