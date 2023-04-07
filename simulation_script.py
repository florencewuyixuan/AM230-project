%matplotlib widget
import gpumd as md
import time

Nsteps = 1000

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
    
    mask = compute_mask(filename)
    field = diffractio.propagate(mask)
    velocity = convertTovVelocity(field)

    update_json(filename, velocity)


