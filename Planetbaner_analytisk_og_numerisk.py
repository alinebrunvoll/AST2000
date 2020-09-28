#  EGEN KODE

'''
Calculating our planet trajectories analytical and numerical
and actually getting a satisfying awnser!!
'''
import numpy as np
import matplotlib.pyplot as plt
import time
import ast2000tools.constants as const
import ast2000tools.utils as utils
from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
from numba import njit
seed = utils.get_seed('alinerb')
system = SolarSystem(seed)
mission = SpaceMission(seed)



def simulate_orbits(system, rot, n_time_steps_per_year):
    G = 4*np.pi**2                      # in Astronomical uits (AU)
    a = system.semi_major_axes
    e = system.eccentricities
    f = np.linspace(0, 2*np.pi, n_time_steps_per_year)
    for i in range(system.number_of_planets):
        r_ = (a[i]*(1-e[i]**2)) / (1 + e[i]*np.cos(np.pi - system.aphelion_angles[i] + f))
        x = r_*np.cos(f)
        y = r_*np.sin(f)
        plt.plot(x, y)
    plt.plot(0, 0, 'yo', label = 'Stjerna')
    plt.legend()
    plt.title('Wolfram Alpha')

    planet_mass = system.masses
    P = np.sqrt(4*np.pi**2/(G*(planet_mass[0]+system.star_mass)) *a[0]**3)
    init_x, init_y = system.initial_positions
    init_vel_x, init_vel_y = system.initial_velocities

    r = np.zeros((n_time_steps_per_year*rot, 7, 2))
    v = np.zeros((n_time_steps_per_year*rot, 7, 2))
    acceleration = np.zeros((n_time_steps_per_year*rot, 7, 2))
    times = np.linspace(0, rot*P, rot*n_time_steps_per_year)
    dt = times[1] - times[0]

    for i in range(system.number_of_planets):
        v_0 = np.array([init_vel_x[i], init_vel_y[i]])
        r_0 = np.array([init_x[i], init_y[i]])
        acceleration[0, i] = -(G*(system.star_mass) / np.linalg.norm(r_0)**3) * r_0
        r[0, i] = r_0
        v[0, i] = v_0

        for t in range(1, n_time_steps_per_year*rot):
            acceleration[t, i] = -(G*(system.star_mass)) /np.linalg.norm(r[t-1, i])**3 * r[t-1, i]
            r[t, i] = r[t-1, i] + v[t-1, i]*dt + 0.5*acceleration[t-1, i]*dt**2
            acceleration[t, i] = -(G*(system.star_mass)) /np.linalg.norm(r[t, i])**3 * r[t, i]
            v[t, i] = v[t-1, i] + 0.5*(acceleration[t-1, i] + acceleration[t, i])*dt


    for planets in range(system.number_of_planets):
        plt.plot(r[:, planets, 0], r[:, planets, 1])
    plt.axis('equal')
    plt.show()

    return times, np.transpose(r), P, dt

# times, r, P, dt = simulate_orbits(system, 25, 10000)



def area(r, P, dt):
    for i in range(system.number_of_planets):
        a1 = 1/2 * r[:, i, 0:1000] * (2*np.pi*r[:, i, 0:1000])/P * dt
        a2 = 1/2 * r[:, i, 5000:6000] * (2*np.pi*r[:, i, 5000:6000])/P * dt
        print(f'Differansen mellom de to arealene til planet {i}: {np.abs(np.sum(a1)-np.sum(a2)):16.14}   |  Den relative feilen: {np.abs(np.sum(a1)-np.sum(a2))/np.abs(np.sum(a2))}')

# area(r, P, dt)

print(time.process_time())

from ast2000tools.shortcuts import SpaceMissionShortcuts
from ast2000tools.space_mission import SpaceMission
# Koden er 74952.
mission = SpaceMission(seed)
shortcuts = SpaceMissionShortcuts(mission, [74952])
# set_launch_parameters() missing 6 required positional arguments: 'thrust', 'mass_loss_rate', 'initial_fuel_mass', 'estimated_launch_duration', 'launch_position', and 'time_of_launch'
shortcuts.thrust
mission.set_launch_parameters(shortcuts.thru)
mission.launch_rocket()
consumed_fuel_mass, final_time, final_position, final_velocity = shortcuts.get_launch_results()
mission.verify_launch_result(final_position)
