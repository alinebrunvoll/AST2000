#  EGEN KODE
import numpy as np
import matplotlib.pyplot as plt

import ast2000tools.constants as const
import ast2000tools.utils as utils
from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
seed = utils.get_seed('alinerb')
system = SolarSystem(seed)

def simulate_orbits(system, duration, n_time_steps_per_year):
    G = 4*np.pi**2            # in Astronomical uits (AU)
    '''
    Here you can implement the orbital simulation for challenge A
    of Part 2 of the project.
    '''
    a = system.semi_major_axes
    e = system.eccentricities
    f = np.linspace(0, 2*np.pi, n_time_steps_per_year)
    for i in range(system.number_of_planets):
        r_ = (a[i]*(1-e[i]**2)) / (1 + e[i]*np.cos(f))
        x = r_*np.cos(f)
        y = r_*np.sin(f)
        plt.plot(x, y)
    plt.plot(0, 0, 'yo')
    plt.title('Our solar system')
    plt.show()


    planet_mass = system.masses
    init_x, init_y = system.initial_positions
    init_vel_x, init_vel_y = system.initial_velocities
    # init_f = system.initial_orbital_angles
    dt = duration/n_time_steps_per_year
    r = np.zeros([n_time_steps_per_year, 7, 2])
    # v_0 = np.zeros([7, 2])
    v = np.zeros([n_time_steps_per_year, 7, 2])
    acceleration = np.zeros([n_time_steps_per_year, 7, 2])
    acceleration_1 = np.zeros([n_time_steps_per_year, 7, 2])
    times = np.zeros([n_time_steps_per_year])

    for t in range(n_time_steps_per_year):
        for i in range(system.number_of_planets):
            v[0, i] = [init_vel_x[i], init_vel_y[i]]
            r[0, i] = [init_x[i], init_y[i]]
            acceleration[t, i, :] = -(G*(system.star_mass + planet_mass[i])) /np.linalg.norm(r[t, i, :])**3 * r[t, i, :]

            # times += dt
            # r[t, i, :] += v[t, i, :]*dt + 0.5*acceleration[t, i, :]*dt**2
            # acceleration_1[t, i, :] = -(G*(system.star_mass + planet_mass[i])) /np.linalg.norm(r[t, i, :])**3 * r[t, i, :]
            # v[t, i, :] += 0.5*(acceleration[t, i, :] + acceleration_1[t, i, :])*dt
            # acceleration[t, i, :] = acceleration_1[t, i, :]

            acceleration[t+1, i, :] = -(G*(system.star_mass + planet_mass[i])) /np.linalg.norm(r[t, i, :])**3 * r[t, i, :]
            r[t+1, i, :] = r[t, i, :] + v[t, i, :]*dt + 0.5*acceleration[t, i, :]*dt**2
            v[t+1, i, :] = v[t, i, :] + 0.5*(acceleration[t, i, :] + acceleration[t+1, i, :])*dt

    plt.plot(r[:, :, 0], r[:, :, 1], 'ro')
    plt.title('Våre planetbaner')
    plt.show()
    # P = np.sqrt(a**3)
    # b = a*np.sqrt(1-e**2)
    # Circ = 2*np.pi * np.sqrt(0.5*(a**2 + b**2))

    # You will probably also need these quantities:
    # const.G_sol
    # system.number_of_planets
    # system.initial_positions
    # system.initial_velocities
    # system.star_mass

    return times, r
simulate_orbits(system,20,1000)
'''
# Prevent the following code from executing when calling `import part_2`
if __name__ == '__main__':

    # Print a message if a newer version of ast2000tools is available
    utils.check_for_newer_version()

    # Construct SpaceMission instance for my mission
    seed = utils.get_seed('alinerb')
    mission = SpaceMission(seed)

    # Extract associated SolarSystem object
    system = mission.system

    # Simulate orbits
    times, planet_positions = simulate_orbits(system,
                                              duration,
                                              n_time_steps_per_year)

    # Verify simulated orbits and write exact orbits to planet_trajectories.npy
    mission.verify_planet_positions(duration,
                                    planet_positions,
                                    filename='planet_trajectories.npy')

    # Generate orbit video to view in SSView
    mission.generate_orbit_video(times,
                                 planet_positions,
                                 filename='orbit_video.xml')
'''
