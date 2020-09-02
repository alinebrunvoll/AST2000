# EGEN KODE

import numpy as np
from scipy import integrate
from sympy import Derivative
from numba import njit
import ast2000tools.utils as utils
import time
# seed = utils.get_seed('XXX')
# 99586
k = 1.38 * 10**(-23)    # Boltzmann-konstanten
m = 2*1.66*10**(-27)    # Mass of hydrogen molecule (kg)
M = 1100

@njit
def combustionchamber(L, T, N, iterations):
    part_pos = np.zeros((iterations, N, 3))
    part_vel = np.zeros((iterations, N, 3))

    init_pos = np.random.uniform(-L/2, L/2, size=(N, 3)) # initial positions randomly generated
    init_vel = np.random.normal(0, np.sqrt(k*T/m), size=(N, 3)) # initial velocities randomly generated

    # Integrating using the Euler-Chromer method

    def motion():
        dt = 10**(-12)
        part_pos[0] = init_pos
        part_vel[0] = init_vel


        for i in range(iterations):
            for j in range(N):
                for k in range(3):
                    if part_pos[i][j][k] <= 0 or part_pos[i][j][k] >= L :
                        part_vel[i][j][k] = - part_vel[i][j][k]


        # Integrating the array with time elements, where each element has N positional vectors in 3 dimentions.
            part_pos[i] = part_pos[i-1] + dt*part_vel[0]
        p_part = 0
        count = 0
        thrust_force = 0
        for i in range(iterations):
            for j in range(N):
                for k in range(3):
                     if part_pos[i, j, 2] <= 0 and np.abs(part_pos[i, j, 0]) < 0.25*L and np.abs(part_pos[i, j, 1]) < 0.25*L:

                        p_part += m*np.abs(part_vel[i, j, 2])
                        count += 1
                        part_pos[i, j, 2] = L


        return part_pos, p_part, count
    return motion()




L = 10**(-6)        # Length of box (m)
T = 3*10**3         # Temperature of gas (K)
N = 10**5           # amt of particles
iterations = 1000   # iterasjons

part_pos, p_part, count = combustionchamber(L, T, N, iterations)

fuel_consumption = count * m   # kg
thrust_force = p_part / 10**(-9)


def fuelconsumed(thrust_force, fuel_consumption, m_rocket, fuel_mass, amt_boxes, delta_v, time, iterations):
    dt = time/iterations






fuelconsumed(thrust_force, fuel_consumption, 1100, 100, xxx, )
