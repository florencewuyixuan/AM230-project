import numpy as np
from Python.pymd.builder import *
from Python.pymd.md import *               # Import the md module from the pymd package
from diffractio import sp, nm, plt, np, mm, degrees, um
from diffractio.scalar_fields_XY import Scalar_field_XY
from diffractio.scalar_sources_XY import Scalar_source_XY
from diffractio.scalar_masks_XY import Scalar_mask_XY
from diffractio.utils_drawing import draw_several_fields

import os
import pdb

from scipy.io import savemat


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
        for ii in range(0,len(particles_x)):
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



def json_to_mat(filename):
     with open(filename) as f:
            data = json.load(f)
            print(filename)
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
                # theta = uniform(-pi,pi)
                # nx, ny = cos(theta), sin(theta)
                # vx, vy = 0.0, 0.0
                # fx, fy = 0.0, 0.0
                # if 'n' in p:  nx, ny = p['n']
                # if 'v' in p:  vx, vy = p['v']
                # if 'f' in p:  fx, fy = p['f']
                particles_x.append(x)
                particles_y.append(y)

            mdic = {"x": particles_x,
                    "y": particles_y}
            # pdb.set_trace()
            os.path.split(filename)
            savemat(os.path.splitext(filename)[0]+'.mat',mdic)

