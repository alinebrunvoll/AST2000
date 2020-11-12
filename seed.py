# IN-BETWEEN PROGRAM WITH LOTS OF INFO FOR ME YES
import numpy as np


# IMPORTS FROM AST2000TOOLS:
import ast2000tools.utils as utils
from ast2000tools.solar_system import SolarSystem
import ast2000tools.constants as const
from ast2000tools.space_mission import SpaceMission

# MAKING THE SEED AND SYSTEM:
seed = utils.get_seed('alinerb')
system = SolarSystem(seed)
mission = SpaceMission(seed)

print(system.print_info())

# VANSKELIG:
'''
import ast2000tools.space_mission as SpaceMission
SpaceMission.get_sky_image_pixel
'''

# INFO ABOUT STAR:
'''
print(f'\n Radius of star: {system.star_radius} km')  # km
print(f'\n Temperature of star: {system.star_temperature} K\n') #K
print(f'\n Mass of star: {system.star_mass} Solar Masses')  # m_sol
'''

# INFO ABOUT PLANETS:
'''
print(system.print_info())
for planet_idx in range(system.number_of_planets):
    print('Planet {:d} is a {} planet with a semi-major axis of {:g} AU.'
         .format(planet_idx, system.types[planet_idx], system.semi_major_axes[planet_idx]))
'''

'''
# SOLAR PANELS AT FRODERIA:
#                  W/m^2/K^4                K                           km -> m                             AU -> m
A = (40/0.12) / (const.sigma * (system.star_temperature)**4 * ( (system.star_radius*1000) / (system.semi_major_axes[3]*const.AU) )**2 )
print(f'\n Area of solar panels on Froderia must be greater than or equal to: {A} m^2')

# SOLAR PANELS ON EARTH:
#                     m             K               m              m
A = (40/0.12) / (const.sigma * (5778)**4 * ( (const.R_sun) / (const.AU) )**2 )
print(f'\n Area of solar panels on Earth must be greater than or equal to: {A} m^2')
'''

'''
# SMALLEST POSSIBLE DISTANCE TO TAKE PICTURE:
P = 2.8*10**(-8)    # km
F = 0.25            # deg
R = system.radii[3] # km

L = (R*P)/F

print(f'must be less than approximately {L} km away from planet.')
'''


# TULL OG TØYS:
'''
x = np.array([4, 5, 6, 7, 8])
y = np.array([5, 2, 3, 1, 4])


plt.plot(x, y)
plt.show()
'''

'''
x = np.matrix([[3, 2], [4, 5]])
print( x * ([5,5]) )
'''

'''
x = np.array([[3, 2], [4, 5]])

print( np.dot(x, (5, 5)) )
'''

'''
list = []

lost = [2,3,4,5,6]

list.append(lost)

last = [7,7,7,7,7]
list.append(last)
print(list)
'''

x = np.zeros(3)
x[:2] = (1,2)

#print(x)
