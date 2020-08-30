# EGEN KODE
import ast2000tools.utils as utils

seed = utils.get_seed('alinerb')

from ast2000tools.solar_system import SolarSystem

system = SolarSystem(seed)



import numpy as np

'''
def f(Thrust, Fuel_Consumption, Init_Mass, Speed_Boost):

    return Consumed_Fuel

Init_Mass = ast2000tools.something()
Thrust =
Fuel_Consumption =
Speed_Boost =
'''

from ast2000tools.space_mission import SpaceMission
mission = SpaceMission(seed)

print(mission.spacecraft_mass)
