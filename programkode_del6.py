# Egen kode
'''
Noe om hva koden gjør
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
import time
from programkode_del4_klasse import did_you_just_assume_my_orientation, launching_sequence
seed = utils.get_seed('alinerb')
mission = SpaceMission(seed)
launch = launching_sequence(mission)
orient = did_you_just_assume_my_orientation(mission)
k = const.k_B
c = 299792458
m_H = const.m_H2/2


timesteps = np.load('times.npy')
planet_positions = np.load('planet_positions.npy')
pos_planet = interpolate.interp1d(timesteps, planet_positions, axis = 0, bounds_error = False, fill_value = "extrapolate")

def read_the_files():
    spectrum_file = open('spectrum_seed64_600nm_3000nm.txt', 'r')
    sigma_noise_file = open('sigma_noise.txt', 'r')

    all_wavelenghts = []
    F_obs = []
    for line in spectrum_file:
        s = line.split()
        all_wavelenghts.append(float(s[0]))
        F_obs.append(float(s[1]))

    sigma_noise = []
    for line in sigma_noise_file:
        s = line.split()
        sigma_noise.append(float(s[1]))

    spectrum_file.close()
    sigma_noise_file.close()
    all_wavelenghts = [i*1e-9 for i in all_wavelenghts]
    return all_wavelenghts, F_obs, sigma_noise

def finding_area_of_interest(all_wavelenghts, F_obs, sigma_noise):
    c = 299792458                      # m/s
    lambda_ideal = [632e-9, 690e-9, 760e-9, 720e-9, 820e-9, 940e-9, 1400e-9, 1600e-9, 1660e-9, 2200e-9, 2340e-9, 2870e-9]
    delta_lambda_max = []
    [delta_lambda_max.append(10000/c * lambda_ideal[i]) for i in range(len(lambda_ideal))]

    areas_of_interest = []
    F_of_interest = []
    sigma_of_interest = []
    for i in range(len(lambda_ideal)):
        area = []
        F_in_area = []
        sigma_in_area = []
        for l in range(len(all_wavelenghts)):
            if (all_wavelenghts[l] > (lambda_ideal[i] - delta_lambda_max[i])) and (all_wavelenghts[l] < (lambda_ideal[i] + delta_lambda_max[i])):
                area.append(all_wavelenghts[l])
                F_in_area.append(F_obs[l])
                sigma_in_area.append(sigma_noise[l])
        areas_of_interest.append(area)
        F_of_interest.append(F_in_area)
        sigma_of_interest.append(sigma_in_area)

    return areas_of_interest, F_of_interest, sigma_of_interest, delta_lambda_max, lambda_ideal

def model_func(areas_of_interest_i, lambda_0, sigma, F_min):
    return (1 + (F_min - 1)*np.exp(-0.5*((areas_of_interest_i - lambda_0)/sigma)**2))

# @njit
def Xi_squared(delta_lambda_max_i, areas_of_interest_i, m_part_i, lambda_ideal_i, F_of_interest_i, sigma_of_interest_i):
    c = 299792458     # m/s
    k = const.k_B
    lambda_0 = np.linspace(lambda_ideal_i - delta_lambda_max_i, lambda_ideal_i + delta_lambda_max_i, 20)
    sigma = np.linspace(lambda_ideal_i/c * np.sqrt(k*150/m_part_i), lambda_ideal_i/c * np.sqrt(k*450/m_part_i), 20)
    F_min = np.linspace(0.7, 1, 20)
    best = 1e100

    for i in range(len(lambda_0)):
        for j in range(len(sigma)):
            for k in range(len(F_min)):
                model = model_func(areas_of_interest_i, lambda_0[i], sigma[j], F_min[k])

                Xi = np.sum(((F_of_interest_i - model)/sigma_of_interest_i)**2)
                if Xi < best:
                    best = Xi
                    lambda_0_best = lambda_0[i]
                    sigma_best = sigma[j]
                    F_min_best = F_min[k]

    return lambda_0_best, sigma_best, F_min_best


all_wavelenghts, F_obs, sigma_noise = read_the_files()
areas_of_interest, F_of_interest, sigma_of_interest, delta_lambda_max, lambda_ideal = finding_area_of_interest(all_wavelenghts, F_obs, sigma_noise)

m_part = [2.6566962e-26, 2.6566962e-26, 2.6566962e-26, 2.988e-26, 2.988e-26, 2.988e-26, 7.3065e-26, 7.3065e-26, 2.66e-26, 2.66e-26, 9.7860872e-26, 7.30637e-26]           # kg
name = ['Oksygen ved bølgelengden 632nm', 'Oksygen ved bølgelengden 690nm', 'Oksygen ved bølgelengden 760nm', 'Vann ved bølgelengden 720nm', 'Vann ved bølgelengden 820nm', 'Vann', 'Karbondioksid', 'Karbondioksid', 'CH4', 'CH4', 'CO', 'N2O']

mean_molecular_weight = (m_part[0] + m_part[3] + m_part[6] + m_part[10])/(4*m_H)


# for i in range(len(areas_of_interest)):
#     lambda_0_best, sigma_best, F_min_best = Xi_squared(delta_lambda_max[i], areas_of_interest[i], m_part[i], lambda_ideal[i], F_of_interest[i], sigma_of_interest[i])
#     T = m_part[i] * sigma_best**2 * c**2 /(k*lambda_ideal[i]**2)
#     v_r = lambda_ideal[i]*c/lambda_0_best # ehhhh nooo
#     print(v_r)
#
#     plt.plot(areas_of_interest[i], F_of_interest[i])
#     # plt.plot(areas_of_interest[i], model_func(areas_of_interest[i], lambda_0_best, sigma_best, F_min_best))
#     plt.title(name[i])
#     plt.xlabel('Bølgelengde [nm]')
#     plt.ylabel('Fluks')
#     plt.show()



# Froderias overflatetemperatur = 215K!

def atmosphere_model(heights, mean_molecular_weight, number_of_steps):
    gamma = 1.4
    G = const.G
    k = const.k_B
    m_H = const.m_H2/2

    # Using the code from part 3
    distance = mission.system.semi_major_axes[3]*const.AU/1000
    T_0 = mission.system.star_temperature*(mission.system.star_radius**2/(4*distance**2))**(1/4)

    T = np.zeros([len(heights)])
    rho = np.zeros([len(heights)])
    P = np.zeros([len(heights)])

    T[0] = T_0
    rho[0] = mission.system.atmospheric_densities[3]
    M = mission.system.masses[3] * const.m_sun
    P[0] = rho[0]*k*T[0]/(mean_molecular_weight*m_H)

    delta_r = (heights[-1]-mission.system.radii[3])/number_of_steps
    C = P[0]**(1-gamma) * T[0]**gamma

    for i in range(len(heights)-1):

        g = G*M/heights[i]**2
        P[i+1] = P[i] + (-rho[i]*g)*delta_r

        if T[i] >= T_0/2:
            T[i+1] = C**(1/gamma) * P[i+1]**(1-(1/gamma))
        else:
            T[i+1] = T_0/2

        rho[i+1] = (P[i+1]*mean_molecular_weight*m_H)/(k*T[i+1])


    return rho, T
print(mission.system.radii[3])
print(2*np.pi/mission.system.rotational_periods[3]*24*60)
number_of_steps = 1000000
heights = np.linspace(mission.system.radii[3]*1000, 8000*1000, number_of_steps)
rho, T = atmosphere_model(heights, mean_molecular_weight, number_of_steps)
print(T)
print(rho)

plt.plot(T, heights, label = 'Temperatur')
plt.legend()
plt.show()

plt.plot(rho, heights, label = 'Atmosfæretetthet')
plt.legend()
plt.show()






print(f'The time to run: {time.process_time()} s')
