# EGEN KODE

import numpy as np
from scipy import integrate

L = 10**(-6) # Length of box (m)
T = 3*10**3 # Temperature of gas (K)
N = 100 # amt of particles



part_pos = np.zeros(N)
part_vel = np.zeros(N)

init_pos = np.random.uniform(0, L+1, size=(N, 3))
init_vel = np.random.uniform(0, L+1, size=(N, 3))





#init_pos = np.array([init_pos_x, init_pos_y, init_pos_z]) # initial positions

'''
for i in range(N+1):
    init_pos = np.array([ np.random.uniform(0, L+1, N)[i], np.random.uniform(0, L+1, N)[i], np.random.uniform(0, L+1, N)[i] ])
    init_pos.append(init_pos[i])

print(init_pos)
'''

'''
def motion():
    dt = 0.001
    v = np.zeros(N)
    x = np.zeros(N)
    t = np.zeros(N)
    #for i in range(N-1):

if (part_pos == 0 or part_pos == L): # Particles hit the walls
'''
