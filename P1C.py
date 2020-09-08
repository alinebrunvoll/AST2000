# EGEN KODE

'''
This code creates a box, and then adds many boxes to create a rocket engine
alternatively it does nothing
'''

import numpy as np
from scipy import integrate
from sympy import Derivative
from numba import njit
import ast2000tools.utils as utils
import time
# seed = utils.get_seed('XXX')
# 99586

k = 1.38 * 10**(-23)                                      # Boltzmann-konstanten
m = 2*1.66*10**(-27)                                      # Mass of hydrogen molecule (kg)
M = 1100

@njit
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
        for j in range(N):
            if part_pos[i, j, 2] <= 0 and np.abs(part_pos[i, j, 0]) < 0.25*L and np.abs(part_pos[i, j, 1]) < 0.25*L:
                p_part += m*np.abs(part_vel[i, j, 2])
                count += 1
                part_pos[i, j, 2] = L

            for xyz in range(3):
                if part_pos[i][j][xyz] <= -L/2 or part_pos[i][j][xyz] >= L/2 :
                    part_vel[i, j, xyz] = - part_vel[i, j, xyz]

                else:
                    # part_vel[i, j, xyz] = (part_pos[i, j, xyz] - part_pos[i-1, j-1, xyz-1])/dt
                    part_vel[i] = (part_pos[i] - part_pos[i-1])/dt
        part_pos[i] +=  part_vel[i]*dt

    return p_part, count



@njit
def fuelconsumed(thrust_force, fuel_consumption, fuel_mass, delta_v):
    mass = M + fuel_mass
    a = thrust_force/(M+fuel_mass)
    delta_t = delta_v / a
    fuel_consumed = fuel_consumption * delta_t
    return fuel_consumed



def rocket_launch(T, L, N, fuel_mass, amt_boxes, dt):
    v_esc = 16592.04                                                  # unnslippningshastighet (m/s)
    g = 14.76                                                         # tyngdeakselerasjonen   (m/s^2)
    p_part, count = combustionchamber(L, T, N, 1000)
    fuel_consumption = (count/10**(-9)) * m * amt_boxes               # kg
    thrust_force = (p_part / 10**(-9) ) * amt_boxes                   # N
    v = [0]
    r = [0]
    t = [0]
    mass = [M + fuel_mass]

    if thrust_force < mass[0]*g:
        print('Thrust force too small!')
        exit()

    i = 0
    while v[i] < v_esc and mass[i] > 1100:
        mass.append(mass[i] - fuel_consumption*dt)
        G = mass[i] * g
        a = (thrust_force - G) / mass[i]
        v.append(v[i] + a*dt)
        r.append(r[i] + v[i+1]*dt)
        t.append(t[i] + dt)
        i += 1

    return print(f'masse:{mass[-1]:.8}, tid (minutter): {t[-1]/60:.5}, posisjon: {r[-1]:.8}, fart: {v[-1]:.8}')


L = 10**(-6)        # Length of box (m)
T = 3*10**3         # Temperature of gas (K)
N = 10**5           # amt of particles
print(combustionchamber(L, T, N, 1000))
# print(rocket_launch(T, L, N, 25*1000, 4*10**15,  1))
print(f'Process time : {time.process_time()}')
