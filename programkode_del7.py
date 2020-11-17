# Egen kode
'''
We are going on a trip in our favorite rocket ship
'''
import numpy as np
import matplotlib.pyplot as plt
import ast2000tools.constants as const
import scipy.interpolate as interpolate
import ast2000tools.utils as utils
from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
from ast2000tools.shortcuts import SpaceMissionShortcuts
from programkode_del4_klasse import did_you_just_assume_my_orientation, launching_sequence
seed = utils.get_seed('alinerb')
mission = SpaceMission(seed)
launch = launching_sequence(mission)
orient = did_you_just_assume_my_orientation(mission)

# Disse lagret vi fra forrige del!
rho_array = np.load('rho.npy')
heights_array = np.load('heights.npy')

rho = interpolate.interp1d(heights_array, rho_array, axis = 0, bounds_error = False, fill_value = "extrapolate")



def air_resistance(A, heights, v):
    rot_period = 2*np.pi*mission.system.radii[3]*1000/(mission.system.rotational_periods[3]*24*60*60)
    C_d = 1

    w = np.zeros([len(heights)])
    F_d = np.zeros([len(heights)])

    for i in range(len(heights)):
        w[i] = rot_period*heights[i]
        v_drag = v - w

        F_d = 1/2 * rho * C_d * A *v_drag**2

        # BURDE VI INTERPOLERE??
    return F_d

def simulating_landing(init_time, init_pos, init_vel, simulation_time, dt):

    G = const.G
    M = mission.system.masses[3] * const.m_sun
    m = mission.lander_mass
    rho_0 = mission.system.atmospheric_densities[3]
    rot_period = 2*np.pi*mission.system.radii[3]*1000/(mission.system.rotational_periods[3]*24*60*60)
    A = mission.lander_area

    N = int(simulation_time/dt)
    vel = np.zeros([N, 2])
    pos = np.zeros([N, 2])
    t = np.zeros([N])

    vel[0] = init_vel
    pos[0] = init_pos
    t[0] = init_time

    for i in range(N-1):

        e_r = pos[i]/np.linalg.norm(pos[i])
        e_theta =  np.array([-e_r[1], e_r[0]])

        g = - G*M/np.linalg.norm(pos[i])**2
        w = rot_period*np.linalg.norm(pos[i]) * e_theta
        v_drag = vel[i] - w

        # F_L = 1/2 * rho_0 * A *(v_t**2 - v_safe**2)    if test og greier
        F_d = 1/2 * rho(np.linalg.norm(pos[i])) * A * np.linalg.norm(v_drag)**2 * (-vel[i]/np.linalg.norm(vel[i]))
        F_g = m*g * e_r
        F_tot = F_d + F_g

        a = F_tot/m

        vel[i+1] = vel[i] + a*dt
        pos[i+1] = pos[i] + vel[i+1]*dt
        t[i+1] = t[i] + dt

        if np.linalg.norm(pos[i]) < mission.system.radii[3]*1000:
            print('Oh no! You crashed!')
            exit()


    return pos, vel, t

# These are the kartesian coordinates printet out form another program
init_time = 100000
init_pos = (-1102189.44814319, -4155404.97443974)
init_vel = (-3825.8392273, 1014.78067287)

simulation_time = 60*60*100
dt = 1

# pos_rocket, vel_rocket, t_rocket = simulating_landing(init_time, init_pos, init_vel, simulation_time, dt)

e_r = init_pos/np.linalg.norm(init_pos)
e_theta =  np.array([-e_r[1], e_r[0]])

init_vel = (-3825.8392273, 1014.78067287)+(e_theta*(100))

pos, vel, t = simulating_landing(init_time, init_pos, init_vel, simulation_time, dt)

# simulation_time = 60*60
# dt = 0.001
# pos_, vel_, t_ = simulating_landing(t[-1], pos[-1], vel[-1], simulation_time, dt)

theta = np.linspace(0, 2*np.pi, 100)
radius = mission.system.radii[3]*1000
figure, axes = plt.subplots(1)
axes.plot(radius*np.cos(theta), radius*np.sin(theta))

plt.plot(pos[:, 0], pos[:, 1])
# plt.plot(pos_[:, 0], pos_[:, 1])
# plt.plot(pos_rocket[:, 0], pos_rocket[:, 1])
plt.axis('equal')
plt.show()

v_terminal = np.sqrt(2*const.G*90*mission.system.masses[3]*const.m_sun/((mission.system.radii[3]*1000)**2*mission.system.atmospheric_densities[3]* mission.lander_area))
print(f'Terminal velocity: {v_terminal} m/s')

area = 2*const.G*90*mission.system.masses[3]*const.m_sun/((mission.system.radii[3]*1000)**2*mission.system.atmospheric_densities[3]* 3**2)
print(f'The area we need the parachute to be: {area} m³')
