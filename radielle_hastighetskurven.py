#  EGEN KODE

'''
Noe om hva funksjonen gjør...
'''
import numpy as np
import matplotlib.pyplot as plt
import ast2000tools.constants as const
import ast2000tools.utils as utils
from ast2000tools.space_mission import SpaceMission
from ast2000tools.solar_system import SolarSystem
seed = utils.get_seed('alinerb')
system = SolarSystem(seed)
