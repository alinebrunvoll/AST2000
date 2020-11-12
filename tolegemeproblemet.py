#  EGEN KODE

'''
Plotting the movement of our star, Universitas, and our home planet, Mørkerius,
(Which turns out to be a straight line...)
'''
import numpy as np
import matplotlib.pyplot as plt
import ast2000tools.constants as const
import ast2000tools.utils as utils
from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
seed = utils.get_seed('alinerb')
system = SolarSystem(seed)

def twobodyproblem(duration, iterations):
    G = 4*np.pi**2                             # in Astronomical units (AU)
    a = system.semi_major_axes
    e = system.eccentricities
    planet_mass = system.masses
    P = np.sqrt((const.G_sol*a[0]**3)/(G*(planet_mass[0]*system.star_mass)))*0.00273973
    init_x, init_y = system.initial_positions
    init_vel_x, init_vel_y = system.initial_velocities

    dt = duration/iterations

    u_pos = np.zeros((iterations, 2))          # Lager initialbetingelser for
    m_pos = np.zeros((iterations, 2))          # planeten vår Mørkerius
    u_vel = np.zeros((iterations, 2))          # og stjerna vår Universitas
    m_vel = np.zeros((iterations, 2))
    u_a = np.zeros((iterations, 2))
    m_a = np.zeros((iterations, 2))
    t = np.linspace(0, duration)
    r_0 = np.array((init_x[0], init_y[0]))
    my = (system.star_mass*planet_mass[0]) / (system.star_mass + planet_mass[0])


    m_vel[0] = np.array((init_vel_x[0], init_vel_y[0]))
    u_vel[0] = np.array((0, 0))
    u_pos[0]  = - (my/system.star_mass)* r_0
    m_pos[0] = (my/planet_mass[0]) * r_0
    u_a[0] = -(G*(system.star_mass + planet_mass[0])) / (np.linalg.norm(u_pos[0] - m_pos[0])**3*system.star_mass) * (u_pos[0] - m_pos[0])
    m_a[0] = -(G*(system.star_mass + planet_mass[0])) / (np.linalg.norm(m_pos[0] - u_pos[0])**3*planet_mass[0]) * (m_pos[0]  - u_pos[0])

    for i in range(1, iterations):
        r = u_pos[i] - m_pos[i]

        u_a[i] = -(G*(system.star_mass + planet_mass[0])) / (np.linalg.norm(u_pos[i-1] - m_pos[i-1])**3*system.star_mass) * (u_pos[i-1]-m_pos[i-1])
        m_a[i] = -(G*(system.star_mass + planet_mass[0])) / (np.linalg.norm(m_pos[i-1] - u_pos[i-1])**3*planet_mass[0]) * (m_pos[i-1] - u_pos[i-1])
        u_pos[i] = u_pos[i-1] + u_vel[i-1]*dt + 0.5*u_a[i-1]*dt**2
        m_pos[i] = m_pos[i-1] + m_vel[i-1]*dt + 0.5*m_a[i-1]*dt**2
        u_a[i] = -(G*(system.star_mass + planet_mass[0])) / (np.linalg.norm(u_pos[i] - m_pos[i])**3*system.star_mass) * (u_pos[i]-m_pos[i])
        m_a[i] = -(G*(system.star_mass + planet_mass[0])) / (np.linalg.norm(m_pos[i]- u_pos[i])**3*planet_mass[0]) * (m_pos[i] - u_pos[i])
        u_vel[i] = u_vel[i-1] + 0.5*u_a[i-1] + u_a[i]*dt
        m_vel[i] = m_vel[i-1] + 0.5*m_a[i-1] + m_a[i]*dt

    plt.plot(m_pos[:, 0], m_pos[:, 1],label = 'Mørkerius')
    plt.plot(u_pos[:, 0], u_pos[:, 1], label = 'Universitas')
    plt.axis('equal')
    plt.legend()
    plt.show()

twobodyproblem(100, 10000)
