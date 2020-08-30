# EGEN KODE

import numpy as np
from scipy import integrate

L = 10**(-6)    # Length of box (m)
T = 3*10**3     # Temperature of gas (K)
N = 100         # amt of particles
n = 100

part_pos = np.zeros((n, N, 3))
part_vel = np.zeros((n, N, 3))


init_pos = np.random.uniform(0, L+1, size=(N, 3)) # initial positions randomly generated
init_vel = np.random.uniform(0, L+1, size=(N, 3)) # initial velocities randomly generated


# Integrating using the Euler-Chromer method
def motion():
    dt = 10**(-9)
    t = np.zeros(n)
    part_pos[0] = init_pos
    part_vel[0] = init_vel
    t[0] = 0

# Integrating the array with 100 elements, where each element has N positional vectors in 3 dimentions.
    for i in range(n):
        for j in range(N):
            for k in range (3):
                if part_pos[i][j][k] <= 0 or part_pos[i][j][k] >= L :
                    part_vel[i][j][k] = - part_vel[i][j][k]

                part_pos[i][j][k] = part_pos[i-1][j-1][k-1] + dt * part_vel[i][j][k]
                t[i] = t[i-1] + dt

            if part_pos[:,:,2] >= (1/2)*L - 0.25/2 and part_pos[:,:,2] <= (1/2)*L + 0.25/2

    return part_pos
