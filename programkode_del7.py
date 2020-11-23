# Egen kode
'''
This code lets us down...
'''
import numpy as np
import matplotlib.pyplot as plt
import ast2000tools.constants as const
import scipy.interpolate as interpolate
import ast2000tools.utils as utils
from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
from ast2000tools.shortcuts import SpaceMissionShortcuts
from numba import njit

seed = utils.get_seed('alinerb')
mission = SpaceMission(seed)


# Disse lagret vi fra forrige del!
rho_array = np.load('rho.npy')
heights_array = np.load('heights.npy')
# Interpolating atmosphere density:
rho = interpolate.interp1d(heights_array, rho_array, axis = 0, bounds_error = False, fill_value = "extrapolate")

# Defining constants:
G = const.G
M = mission.system.masses[3] * const.m_sun
m = mission.lander_mass
rho_0 = mission.system.atmospheric_densities[3]
rot_vel = 2*np.pi/(mission.system.rotational_periods[3]*24*60*60)
A = mission.lander_area
R = mission.system.radii[3]*1000
v_safe = 3

# Calculating terminal velocity and area needed for parachute:
v_terminal = np.sqrt( 2 * G * 90 * M/(R**2 * rho_0 * A) )
Chute_area = 2 * G * 90 * M / (R**2 * rho_0 * 3**2)
#print(f'Terminal velocity: {v_terminal} m/s')
#print(f'The area we need the parachute to be: {Chute_area} m³')


#@njit
def simulating_landing(init_time, init_pos, init_vel, simulation_time, dt):

    N = int(simulation_time/dt)
    vel = np.zeros([N, 2])
    pos = np.zeros([N, 2])
    t = np.zeros([N])
    A = mission.lander_area

    vel[0] = init_vel
    pos[0] = init_pos
    t[0] = init_time

    chute = 0; thrust = 0
    thrust_point = 0; chute_point = 0; land_point = 0

    # Integration loop
    for i in range(N-1):

        e_r = pos[i]/np.linalg.norm(pos[i])
        e_theta =  np.array([-e_r[1], e_r[0]])

        g = G*M/np.linalg.norm(pos[i])**2
        w = rot_vel*np.linalg.norm(pos[i]) * e_theta
        v_drag = vel[i] - w

        F_d = 1/2 * rho(np.linalg.norm(pos[i])) * A * np.linalg.norm(v_drag)**2 * (-v_drag/np.linalg.norm(v_drag))
        F_g = - m*g * e_r
        F_L = 1/2 * rho_0 * A *(v_terminal**2 - v_safe**2) * (e_r)

        # Landing Thruster:
        if np.linalg.norm(pos[i]) <= (R + 200):
            F_tot = F_d + F_g + F_L
            if thrust == 0:
                print('Landing thrusters activated\n')
                thrust_point = pos[i]
                thrust = 1
        else:
            F_tot = F_d + F_g

        a = F_tot/m

        vel[i+1] = vel[i] + a*dt
        pos[i+1] = pos[i] + vel[i+1]*dt
        t[i+1] = t[i] + dt

        # Hvis vi krasjer i planeten:
        if np.linalg.norm(pos[i]) <= R:
            print('Touchdown!\n')
            if np.linalg.norm(vel[i]) > 3:
                print('You crashed though...\n')
            land_point = pos[i]

        # Hvis vi brenner opp i atmosfæren:
        if np.linalg.norm(F_d) > (10**7 * A):
            print('Sorry... Your lander just couldn\'t handle the pressure...\n')

        # Fallskjerm:
        if np.linalg.norm(pos[i]) <= R + 100000:
            A = mission.lander_area + Chute_area
            if chute == 0:
                chute_point = pos[i]
                print('\nParachute deployed\n')
                chute = 1
            # Hvis fallskjermen ryker:
            if np.linalg.norm(F_d) > 250000:
                print('Oh, chute! Looks like your parachute was caught by the wind...\n')

    return pos, vel, t, chute_point, thrust_point, land_point


# These are the cartesian coordinates printed out from another program:
init_time = 100000
init_pos = (4127634.61920073, 1202017.45199339)
init_vel = (-1106.69123521,  3800.27100878)

# Unit vectors:
e_r = init_pos/np.linalg.norm(init_pos)
e_theta = np.array([-e_r[1], e_r[0]])

# Simulation:
simulation_time = 60*60*3.18
dt = 0.01
init_vel = init_vel - (e_theta*(75))
pos, vel, t, chute_point, thrust_point, land_point = simulating_landing(init_time, init_pos, init_vel, simulation_time, dt)

# Plotting planet radius:
theta = np.linspace(0, 2*np.pi, 100)
radius = mission.system.radii[3]*1000
figure, axes = plt.subplots(1)
axes.plot(radius*np.cos(theta), radius*np.sin(theta), label='Froderia radius')

# Plotting lander position:
plt.plot(pos[:, 0], pos[:, 1],label='Landingsmodul')
plt.axis('equal')
plt.title('Landing')
plt.xlabel('Avstand [m]')
plt.ylabel('Avstand [m]')
plt.legend()
'''
# Plotting parachute deployment, thruster activation and landing position:
plt.plot(chute_point[0], chute_point[1], 'bH', label='Fallskjerm utløst')
plt.plot(thrust_point[0], thrust_point[1], 'b^', label='Thruster aktivert')
plt.plot(land_point[0], land_point[1], 'r*', label='Touchdown')'''
plt.show()
