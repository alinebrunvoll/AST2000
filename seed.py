import ast2000tools.utils as utils

seed = utils.get_seed('alinerb')

from ast2000tools.solar_system import SolarSystem

system = SolarSystem(seed)


solar_mass = system.star_mass
solar_radius = system.star_radius

print(solar_mass)
print(solar_radius)  # km

for planet_idx in range(system.number_of_planets):
    print('Planet {:d} is a {} planet with a semi-major axis of {:g} AU.'
          .format(planet_idx, system.types[planet_idx], system.semi_major_axes[planet_idx]))



print(system.print_info())
