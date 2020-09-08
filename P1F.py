# EGEN KODE


'''noe om hva koden gjør'''

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

'''
@njit
def combustionchamber(L, T, N, iterations):
    part_pos = np.zeros((iterations, N, 3))
    part_vel = np.zeros((iterations, N, 3))

    init_pos = np.random.uniform(-L/2, L/2, size=(N, 3))          # initial positions randomly generated
    init_vel = np.random.normal(0, np.sqrt(k*T/m), size=(N, 3))   # initial velocities randomly generated

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

        return p_part, count
    return motion()
'''


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
        part_vel[i] = (part_pos[i] - part_pos[i-1])/dt
        part_pos[i] +=  part_vel[i]*dt

        np.where( np.less_equal(part_pos[i,:,2], 0) == True and np.less_equal(part_pos[i, :, 0], 0.25*L) == True, p_part += m*np.abs(part_vel[:, :, 2]) and count += 1 and part_pos[i, :, 2] = L)
            # p_part += m*np.abs(part_vel[:, :, 2])
            # count += 1
            # part_pos[i, :, 2] = L

        for xyz in range(3):
            if np.less_equal(part_pos[i, :, xyz], -L/2) or np.greater_equal(part_pos[i, :, xyz], L/2):
                part_vel[i, :, xyz] *= -1
        # lost_par =
        # new_indices = np.where(lost_par == True)
        #
        # collision_points = np.logical(cond2, cond2)
        # v[collision_points] *= -1


            # if part_pos[:][:][xyz] <= -L/2 or part_pos[:][:][xyz] >= L/2 :
            #     part_vel[:, :, xyz] = - part_vel[:, :, xyz]

    # if part_pos[:, :, 2] <= 0 and np.abs(part_pos[:, :, 0]) < 0.25*L and np.abs(part_pos[:, :, 1]) < 0.25*L:
    #     p_part += m*np.abs(part_vel[:, :, 2])
    #     count += 1
    #     part_pos[:, :, 2] = L


    for i in range(iterations):
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
    v_rot = (2*np.pi/(1.0025*24*60*60)) * 9329.67*1000                # speed from the rotation of the planet (m/s)
    p_part, count = combustionchamber(L, T, N, 1000)
    fuel_consumption = count * m * amt_boxes                          # kg
    thrust_force = (p_part / 10**(-9) ) * amt_boxes                   # N
    v = np.array([[0, v_rot]])
    r = np.array([[0, 0]])
    t = np.array([0])
    mass = np.array([M + fuel_mass])

    # v = [[0, v_rot]]
    # r = [[0, 0]]
    # t = [0]
    # mass = [M + fuel_mass]


    if thrust_force > mass[0]*g:
        print(':)')

    i = 0
    print(v[0][0])
    print(v[0][1])
    while np.sqrt(v[i][0]**2 + v[i][1]**2) < v_esc and mass[i] > 1100:
        mass[i+1] = mass[i] - fuel_consumption*dt
        G = mass[i] * g
        a = [(thrust_force - G) / mass[i], 0]
        v[i] += a*dt
        r[i] += v[i]*dt
        t[i] += dt

        v[i+1] = v[i] + a*dt
        r[i+1] = r[i] + v[i+1]*dt
        t[i+1] = t[i] + dt
        # i += 1

        # v.append(v[i] + a*dt)
        # r.append(r[i] + v[i+1]*dt)
        # t.append(t[i] + dt)
        i += 1



    return print(f'masse:{mass[-1]}, tid (minutter): {t[-1]/60}, posisjon: {r[-1]}, fart: {v[-1]}')


L = 10**(-6)        # Length of box (m)
T = 3*10**3         # Temperature of gas (K)
N = 10**5           # amt of particles

print(combustionchamber(L, T, 20, 5))
# print(rocket_launch(T, L, N, 25*1000, 3*10**15,  1))
print(time.process_time())
