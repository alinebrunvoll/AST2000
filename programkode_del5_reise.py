# Egen kode
'''
We are going on a trip in our favorite rocket ship
'''
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
import ast2000tools.constants as const
import ast2000tools.utils as utils
from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
from ast2000tools.shortcuts import SpaceMissionShortcuts
from programkode_del4_klasse import did_you_just_assume_my_orientation, launching_sequence
seed = utils.get_seed('alinerb')
mission = SpaceMission(seed)
launch = launching_sequence(mission)
orient = did_you_just_assume_my_orientation(mission)
Ms = mission.system.star_mass
M_planet = mission.system.masses
G = const.G_sol

planet_positions = np.load('planet_positions.npy')
planet_velocities = np.load('planet_velocities.npy')
timesteps = np.load('times.npy')
plan_v_rocket = np.load('optimal_rakett_hastighet.npy')
plan_r_rocket = np.load('optimal_rakett_posisjon.npy')
plan_t_rocket = np.load('optimal_rakett_tid.npy')

pos_planet = interpolate.interp1d(timesteps, planet_positions, axis = 0, bounds_error = False, fill_value = "extrapolate")
vel_planet = interpolate.interp1d(timesteps, planet_velocities, axis = 0, bounds_error = False, fill_value = "extrapolate")
ideal_pos = interpolate.interp1d(plan_t_rocket, plan_r_rocket, axis = 0, bounds_error = False, fill_value = "extrapolate")
ideal_vel = interpolate.interp1d(plan_t_rocket, plan_v_rocket, axis = 0, bounds_error = False, fill_value = "extrapolate")

when_to_start = 0.5583448304020219
delta_v1 = 1.6200938839719106
time_of_flight = 0.29420809716565316

'''Launching and cheching if orientation works '''
final_time, final_position, final_velocity = launch.launch_results(when_to_start, np.pi/2)
mission.take_picture(filename='unknown_pic.png')
mission.verify_manual_orientation(orient.where_we_at(final_time, mission.measure_distances()), orient.how_fast_we_goin(mission.measure_star_doppler_shifts()), orient.what_we_looking_at_here())

mission.verbose = False
travel = mission.begin_interplanetary_travel()
current_time, current_pos, current_vel = travel.orient()


'''Making a plotting function'''
def plot_me(t, pos):
    plt.plot(pos_planet(t)[:, 0, 0], pos_planet(t)[:, 0, 1], 'g', label = 'Mørkerius')
    plt.plot(pos_planet(t)[:, 3, 0], pos_planet(t)[:, 3, 1], 'tab:orange', label = 'Froderia')
    plt.plot(pos[0], pos[1], 'r')


# travel.look_in_direction_of_planet(0)
# travel.take_picture(filename='goodbye_home.xml')

'''Adding the delta_v1 from our original plan'''
delta_v1 = current_pos /np.linalg.norm(current_pos) * (delta_v1 + 0.16)
travel.boost(delta_v1)

pos = np.zeros([2,300])
t = np.zeros([300])
for i in range(300):
    travel.coast(0.0001)
    current_time, current_pos, current_vel = travel.orient()
    pos[:, i] = current_pos
    t[i] = current_time
plot_me(t, pos)

# travel.look_in_direction_of_motion()
# travel.take_picture(filename='where_we_at.xml')


'''Adding a boost to correct our cource'''
vel_wrong = current_vel - ideal_vel(current_time)
travel.boost(-vel_wrong)
travel.coast(0.01)

pos = np.zeros([2,1200])
t = np.zeros([1200])
for i in range(1200):
    travel.coast(0.0001)
    current_time, current_pos, current_vel = travel.orient()
    pos[:, i] = current_pos
    t[i] = current_time
plot_me(t, pos)
# plt.plot(plan_r_rocket)
# plt.show()
# sys.exit()

'''
Making a new boost that points towards our planet, and changing it
slightly so that we will go towards where the planet is going to be
'''
# delta_v2 = current_vel /np.linalg.norm(current_vel) * (-0.86) #[0.68544445 0.51938993]
delta_v2 = np.array([0.68544445, 0.303538993])
travel.boost(delta_v2)

# travel.look_in_direction_of_planet(3)
# travel.start_video()
pos = np.zeros([2,432])
t = np.zeros([432])
for i in range(432):
    travel.coast(0.0001)
    current_time, current_pos, current_vel = travel.orient()
    pos[:, i] = current_pos
    t[i] = current_time
# travel.finish_video()
plot_me(t, pos)

# travel.look_in_direction_of_planet(3)
# travel.take_picture(filename='its_Froderia.xml')

'''Finding out if we are close enough to our planet '''
l = np.linalg.norm(current_pos)*np.sqrt(M_planet[3]/(10*Ms))
print(np.linalg.norm(current_pos - pos_planet(t)[-1, 3]) <= l)



'''Making our orbital injection manouever'''
distance = - current_pos + pos_planet(t)[-1, 3]
v_stable = np.sqrt(G*M_planet[3]/np.linalg.norm(distance))
r_hatt = distance/(np.linalg.norm(distance))
theta_hatt = np.array([-r_hatt[1], r_hatt[0]])

v_orbit = - v_stable*theta_hatt + vel_planet(t)[-1, 3]
travel.boost(v_orbit-current_vel)

'''Checking if we still are in orbit '''
pos = np.zeros([2, 200])
t = np.zeros([200])
for i in range(200):
    travel.coast(0.00001)
    current_time, current_pos, current_vel = travel.orient()
    pos[:, i] = current_pos
    t[i] = current_time

plot_me(t, pos)

l = np.linalg.norm(current_pos)*np.sqrt(M_planet[3]/(10*Ms))
print(np.linalg.norm(current_pos - pos_planet(t)[-1, 3]) <= l)


'''Plotting the last point '''
plt.plot(current_pos[0], current_pos[1], 'ro')
plt.plot(pos_planet(t)[-1, 3, 0], pos_planet(t)[-1, 3, 1], 'go', label = 'Froderia')
plt.axis('equal')
plt.xlabel('AU')
plt.ylabel('AU')
plt.show()

travel.record_destination(3)

distance = (np.transpose(current_pos) - pos_planet(t)[-1, 3])


pos = np.zeros([2, 200])
t = np.zeros([200])
for i in range(200):
    travel.coast(0.0001)
    current_time, current_pos, current_vel = travel.orient()
    pos[:, i] = current_pos
    t[i] = current_time
distance = (np.transpose(current_pos) - pos_planet(t)[-1, 3])


'''Finding everything we need to find in task D '''
dist = (np.transpose(pos) - pos_planet(t)[:, 3])

v_theta = dist * np.linalg.norm(current_vel)/np.linalg.norm(dist)
v_r = np.linalg.norm(dist)

my = M_planet[3]*(1100/const.m_sun)/(M_planet[3] + (1100/const.m_sun))
r1_cm = -my/M_planet[3] *distance
r2_cm = my/(1100/const.m_sun) *distance

plt.plot(0, 0, 'bo')
plt.plot(dist[:, 0], dist[:, 1])

periapsis = np.amin(dist)
apoapsis = np.amax(dist)

a = (np.abs(periapsis) + np.abs(apoapsis))/2
P = np.sqrt(4*np.pi**2*a**3/(G*(M_planet[3]+1100/const.m_sun)))
c = apoapsis - a
e = c/a
b = a*np.sqrt(1-e**2)

plt.axis('equal')
plt.xlabel('AU')
plt.ylabel('AU')
plt.show()


''''''''''''''''' ----------- Part 6 ------------------'''''''''''''''''

landing = mission.begin_landing_sequence()
current_time, current_pos, current_vel = landing.orient()

r_hatt = current_pos/(np.linalg.norm(current_pos))
theta_hatt = np.array([-r_hatt[1], r_hatt[0], 0])


landing.boost(theta_hatt *(-350))

'''
# orbital injection manouever
current_time, current_pos, current_vel = landing.orient()
v_stable = np.sqrt(G*M_planet[3]/np.linalg.norm(current_pos))
r_hatt = distance/(np.linalg.norm(current_pos))
theta_hatt = np.array([-r_hatt[1], r_hatt[0], 0])

vel_plan = np.zeros([])
vel_plan[:2] = vel_planet(t)[-1, 3]
v_orbit = - v_stable*theta_hatt + vel_plan
landing.boost(v_orbit-current_vel)
'''

landing.look_in_direction_of_planet(3)
landing.start_video()
pos = np.zeros([20000, 3])
t = np.zeros([20000])
for i in range(20000):
    landing.fall(60*30)
    current_time, current_pos, current_vel = landing.orient()
    pos[i] = current_pos
    t[i] = current_time
'''
    if i == 12183:
        v_stable = np.sqrt(G*M_planet[3]/np.linalg.norm(current_pos))
        r_hatt = current_pos/(np.linalg.norm(current_pos))
        theta_hatt = np.array([-r_hatt[1], r_hatt[0], 0])

        v_orbit = - v_stable*theta_hatt
        print(v_orbit)
        landing.boost(v_orbit-current_vel)
'''

landing.finish_video(filename = 'closeup2.xml')

print(pos)
print(np.shape(pos))

pos_norm = np.zeros([len(pos)])
for i in range(len(pos)):
    pos_norm[i] = np.linalg.norm(pos[i])

print(pos_norm)
print(np.argmin(pos_norm))

print(pos_norm[np.argmin(pos_norm)])
'''
# IKKE SLETT MEG
def coordinates(landing_pos, time_elapsed):
    rot_period = mission.system.rotational_periods[3]
    r = mission.system.radii[3]

    x, y, z = landing_pos
    phi = np.arctan(y/x)
    phi_rot = 2*np.pi*time_elapsed/rot_period

    phi_new = phi + phi_rot

    return (r, phi_new, 0)

current_time, current_pos, current_vel = landing.orient()
landing_pos = current_pos/np.linalg.norm(current_pos) * mission.system.radii[3]

coordinates(landing_pos, current_time)
'''
