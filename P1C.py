# EGEN KODE

import numpy as np
from scipy import integrate
from numba import njit
import ast2000tools.utils as utils
# seed = utils.get_seed('XXX')
# 99586


def combustionchamber(L, T, N, iterations):
    k = 1.38 * 10**(-23)
    m = 1.67*10**(-27)
    part_pos = np.zeros((iterations, N, 3))
    part_vel = np.zeros((iterations, N, 3))

    init_pos = np.random.uniform(-L/2, L/2, size=(N, 3)) # initial positions randomly generated
    init_vel = np.random.normal(0, np.sqrt(k*T/m), size=(N, 3)) # initial velocities randomly generated

    # Integrating using the Euler-Chromer method

    def motion():
        dt = 10**(-12)
        part_pos[0] = init_pos
        part_vel[0] = init_vel

        def particles():
            for i in range(iterations):
                for j in range(N):
                    for k in range(3):
                        if part_pos[i][j][k] <= 0 or part_pos[i][j][k] >= L :
                            part_vel[i][j][k] = - part_vel[i][j][k]
                     # newvel = np.where(np.abs(part_pos[i,j]) >= L/2, -1, 1)
                     # part_vel[i,j] *= newvel
        # Integrating the array with time elements, where each element has N positional vectors in 3 dimentions.
                part_pos[i] = part_pos[i-1] + dt*part_vel[0]
        particles()

        count = 0

        if part_pos[i, j, 2] == 0 and np.abs(part_pos[:, :, 0:1]) < 0.25*L:
            count += 1

        # hole = np.where(part_pos[:, :, 2] = 0 and np.abs(part_pos[:, :, 0]) < 0.25*L and np.abs(part_pos[:, :, 1]) < 0.25*L, True, None)
        # count.append(hole)
        return part_pos

    return motion()


L = 10**(-6)        # Length of box (m)
T = 3*10**3         # Temperature of gas (K)
N = 10**5           # amt of particles
iterations = 10    # time

print(combustionchamber(L, T, N, iterations))
