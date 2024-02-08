# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 14:10:07 2023

@author: JulienJacquemin
"""

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/LaboThapLibrary')

import math
import numpy as np
from scipy.optimize import fsolve, least_squares
import scipy.special as bessel
import random
import csv
import matplotlib.pyplot as plt
from Port.Mass_connector import Mass_connector
from Aerotherme_model import Dry_cooler


"Connectors"
Water_in = Mass_connector()
Water_out = Mass_connector()

Air_in = Mass_connector()
Air_out = Mass_connector()

Water_in.set_fluid('water')
Water_in.set_Q_dot(0.003) # Volume flow rate [m^3/s]
Water_in.set_T(65+273.15) # Water entry temperature

Air_in.set_fluid('air')
Air_in.set_Q_dot(106500/3600/4) # Volume flow rate [m^3/s]
Air_in.set_T(20+273.15) # Air entry temperature
    
"Intput variables"
Fan1 = True # If the fan is working or not
Fan2 = True
Fan3 = True
Fan4 = False

fan = [Fan1, Fan2, Fan3, Fan4]


HX = Dry_cooler()

HX.set_geometrical_parameters(**{
    'nb_row': 25, # Number of row of tubes
    'nb_column': 2, # Number of column of tubes
    'D_ext': 0.0127, # External diameter [m]
    'th_wall': 0.002, # Wall thikness of the tubes [m]
    'e_wall': 0.0001, # Wall roughness [m] (hypothesis) 
    'L_tube': 1.4, # Length of one tube [m]
    'k_tube': 45, # Conductivity of steel [W/mK] (hypothesis on the material of the fins) 

    # "Fins data for one fictional exchanger"
    'th_fin': 0.00045, # Thickness of the fins [m] 
    'p_fin': 0.002217, # space between fins [m] (evaluated by counting the number of fins on 4cm and applying p_fin = (4-19*th_fin)/(nb_fins_4-1))
    'h_fin': 0.17/2, # Height of the fins [m] (hypothesis, because two tubes of different temperature use the same fin, we divided the fin hight by two. 
    # It derives from the hypothesis that the fin is high enough so that the two heat transfers are decoupled)
    'l_fin': 1.2, # Length of the fins [m]
    'k_alu': 237, # Conductivity of aluminium [W/mK] (hypothesis on the material of the fins)
})

HX.set_calibrated_parameters(**{
    # 'C1' : 1, # Coefficient of the heat transfer between the water and the tube
    # 'C2' : 0, # Coefficient of the heat transfer between the air and the fin
    # 'C3' : 1, # Coefficient of the heat transfer between the fin and the tube
    # 'C4' : 0, # Coefficient of the heat transfer between the air and the tube
})

HX.update_inputs(Water_in, Air_in, Water_out, Air_out, fan)

HX.solve()

print(HX.Hxs[1].Tw_out)