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


@njit
def fuelconsumed(thrust_force, fuel_consumption, fuel_mass, delta_v):
    mass = M + fuel_mass
    a = thrust_force/(M+fuel_mass)
    delta_t = delta_v / a
    fuel_consumed = fuel_consumption * delta_t
    return fuel_consumed


# print(fuelconsumed(thrust_force, fuel_consumption, 1100, 1))





L = 10**(-6)        # Length of box (m)
T = 3*10**3         # Temperature of gas (K)
N = 10**5           # amt of particles

@njit
def rocket_launch(T, L, N, fuel_mass, amt_boxes, iterations):
    v_esc = 9097.12           # unnslippningshastighet (m/s)

    part_pos, p_part, count = combustionchamber(L, T, N, iterations)
    fuel_consumption = count * m * amt_boxes                          # kg
    thrust_force = (p_part / 10**(-9) ) * amt_boxes                   # N

    v = np.zeros(iterations)
    r = np.zeros(iterations)
    t = np.zeros(iterations)
    mass = np.zeros(iterations)
    dt = 0.01

    v[0] = 0
    r[0] = 0
    t[0] = 0
    mass[0] = M + fuel_mass


    for i in range(iterations):
        while v[i] < v_esc:
            a = (thrust_force + v[i]*(mass[i+1] - mass[i])) / mass[0]
            v[i] = v[i-1] + a*dt
            r[i] = r[i-1] + v[i]*dt
            t[i] = t[i-1] + dt

print(rocket_launch(T, L, N, 20*1100, 10*12, 500))
