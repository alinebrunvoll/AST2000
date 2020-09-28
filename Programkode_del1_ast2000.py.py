# EGEN KODE
'''
This code creates a box, and then adds many boxes to create a rocket engine!
Alternatively it does nothing...
'''

import numpy as np
from scipy import integrate
from sympy import Derivative
from numba import njit
import ast2000tools.constants as const
import ast2000tools.utils as utils
import time
from ast2000tools.shortcuts import SpaceMissionShortcuts
from ast2000tools.space_mission import SpaceMission
seed = utils.get_seed('alinerb')
mission = SpaceMission(seed)
shortcuts = SpaceMissionShortcuts(mission, [74952])


k = 1.38 * 10**(-23)                                      # Boltzmann-konstanten
# m = 2*1.66*10**(-27)                                      # Mass of hydrogen molecule (kg)
m = const.m_H2
M = 1100


def combustionchamber(L, T, N, iterations):
    part_pos = np.zeros((iterations, N, 3))
    part_vel = np.zeros((iterations, N, 3))
    init_pos = np.random.uniform(-L/2, L/2, size=(N, 3))          # initial positions randomly generated
    init_vel = np.random.normal(0, np.sqrt(k*T/m), size=(N, 3))   # initial velocities randomly generated

    dt = 10**(-12)
    part_pos[0] = init_pos
    part_vel[0] = init_vel
    p_part = 0
    count = 0

    for i in range(iterations):
        part_pos[i] +=  part_vel[i]*dt

# Making the hole
        where_part = np.logical_and(np.less_equal(part_pos[i,:,2], 0), np.logical_and(np.less_equal(np.abs(part_pos[i, :, 0]), 0.25*L) , np.less_equal(np.abs(part_pos[i, :, 1]), 0.25*L)))
# Calculating numper of escaped particles, finding the momentum and putting them back in
        p_part += m*np.sum(part_vel[i, where_part, 2])
        count += np.sum(where_part)
        part_pos[i, where_part, 2] = np.random.uniform(-L/2, L/2)

# Checking for wall colitions
        part_pos_hit_box = np.abs(part_pos[i, :, :]) >= L/2
        part_vel[i, part_pos_hit_box] *= -1
    return -p_part, count

@njit
def fuelconsumed(thrust_force, fuel_consumption, fuel_mass, delta_v):
    mass = M + fuel_mass
    a = thrust_force/(M+fuel_mass)
    delta_t = delta_v / a
    fuel_consumed = fuel_consumption * delta_t
    return fuel_consumed


def rocket_launch(T, L, N, fuel_mass, amt_boxes, dt):
    gamma = 6.67e-11
    R = 9329.67*1000                                           # radius (m)
    M_planet = 1.925e25                                        # kg

    p_part, count = combustionchamber(L, T, N, 1000)
    fuel_consumption = count * m * amt_boxes                   # kg
    thrust_force = (p_part/10e-9) * amt_boxes                  # N

    v = 0
    r = 0
    t = 0
    v_esc = np.array([np.sqrt((2*gamma*M_planet)/R)])           # m/s
    mass = M + fuel_mass

# Checking if thrust_force is bigger than the gravitational pull
    if thrust_force < mass*((gamma*M_planet)/R**2):
        print('Thrust force too small!')
        exit()

# Calculating the rocets Mass, time and speed during liftoff
    i = 0
    while v < v_esc and mass > 1100:
        mass =  mass - fuel_consumption*dt
        g = (gamma*M_planet) / (R+r)**2
        G = mass * g
        a = (thrust_force - G) / mass
        v += a*dt
        r += v*dt
        t += dt
        v_esc = np.sqrt( (2*gamma*M_planet)/(R+r) )
        i += 1

    # return print(f'Masse:{mass}, Tid (minutter): {t/60}, Posisjon: {r}, Fart: {v}')
    return thrust_force, fuel_consumption, fuel_mass, t

L = 10**(-7)        # Length of box (m)
T = 3*10**3         # Temperature of gas (K)
N = 10**5           # amt of particles

thrust_force, fuel_consumption, fuel_mass, t = rocket_launch(T, L, N, 15*1000, 10**18,  1)

# Koden er 74952.
print(shortcuts)

mission.set_launch_parameters(thrust_force, fuel_consumption, fuel_mass, t, (0, 0), 0)
mission.launch_rocket()
consumed_fuel_mass, final_time, final_position, final_velocity = shortcuts.get_launch_results()
mission.verify_launch_result(final_position)
# print(time.process_time())
