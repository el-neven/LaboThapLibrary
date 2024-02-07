# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 10:38:56 2023

@author: JulienJacquemin
"""

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/Python_Library')

from Aerotherme_model import Dry_cooler
import csv
import numpy as np
import random
from scipy.optimize import least_squares
import math
from Port.Mass_connector import Mass_connector

# dry cooler class
HX = Dry_cooler()

HX.set_geometrical_parameters(**{
    'nb_row': 25, # Number of row of tubes
    'nb_column': 2, # Number of column of tubes
    'D_ext': 0.0127, # Wall thikness of the tubes [m]
    'th_wall': 0.002, # Internal diameter [m]
    'D_int': 0.0127-2*0.002, # Internal diameter [m]
    'R_int': 0.0127/2-0.002, # Internal radius [m]
    'e_wall': 0.0001, # Wall roughness [m] (hypothesis) 
    'L_tube': 1.4, # Length of one tube [m]
    'L_tot': 1.4*2*25, # Total cross sectional area of the tubes [m^2]
    'Ac_tot': 25*2*math.pi*(0.0127/2-0.002)**2, # Total cross sectional area of the tubes [m^2]
    'k_tube': 45, # Conductivity of steel [W/mK] (hypothesis on the material of the fins) 

    # "Fins data for one fictional exchanger"
    'th_fin': 0.00045, # Thickness of the fins [m] 
    'p_fin': 0.002217, # space between fins [m] (evaluated by counting the number of fins on 4cm and applying p_fin = (4-19*th_fin)/(nb_fins_4-1))
    'h_fin': 0.17/2, # Height of the fins [m] (hypothesis, because two tubes of different temperature use the same fin, we divided the fin hight by two. 
    # It derives from the hypothesis that the fin is high enough so that the two heat transfers are decoupled)
    'l_fin': 1.2, # Length of the fins [m]
    'k_alu': 237, # Conductivity of aluminium [W/mK] (hypothesis on the material of the fins)
})

"Connectors"
Water_in = Mass_connector()
Water_out = Mass_connector()

Air_in = Mass_connector()
Air_out = Mass_connector()

Water_in.set_fluid('water')
Water_in.set_Q_dot(0.003) # Volume flow rate [m^3/s]
Water_in.set_T(65+273.15) # Water entry temperature

Air_in.set_fluid('air')

def optimize_HX(mode=1):
    if mode == 1:
        prob1 = least_squares(optimize_Hx_tot, [1, 0, 1, 0], bounds=([0, float("-inf"), 0, float("-inf")], [float("inf"), float("inf"), float("inf"), float("inf")]))
        C1 = prob1.x[0]
        C2 = prob1.x[1]
        C3 = prob1.x[2]
        C4 = prob1.x[3]
    if mode == 2:
        prob1 = least_squares(optimize_Hx_ph1, [1, 0], bounds=([0, float("-inf")], [float("inf"), float("inf")]))
        C1 = prob1.x[0]
        C2 = prob1.x[1]
        prob2 = least_squares(optimize_Hx_ph2, [1, 0], args=(C1, C2), bounds=([0, float("-inf")], [float("inf"), float("inf")]))
        C3 = prob2.x[0]
        C4 = prob2.x[1]
    
    return [C1, C2, C3, C4]

def optimize_Hx_tot(x):
    
    C1 = x[0]
    C2 = x[1]
    C3 = x[2]
    C4 = x[3]
    
    HX.set_calibrated_parameters(**{
    'C1' : C1, # Coefficient of the heat transfer between the water and the tube
    'C2' : C2, # Coefficient of the heat transfer between the air and the fin
    'C3' : C3, # Coefficient of the heat transfer between the fin and the tube
    'C4' : C4, # Coefficient of the heat transfer between the air and the tube
    })


    with open("C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/Python_Library/Component/Heat_Exchanger/Aerotherme_model/fake_sample.csv", "r") as fake_file:
        samples = csv.reader(fake_file, delimiter=";")
        diff = []
    
        for i, sample in enumerate(samples):
            if i==0:
                continue
            Tw_out_exp = float(sample[0])
            Tw_in = float(sample[1])
            Ta = float(sample[2])
            Qw = float(sample[3])
            Air_in.set_Q_dot(106500/3600/4) # Volume flow rate [m^3/s]
            Air_in.set_T(Ta) # Air entry temperature

            Water_in.set_Q_dot(Qw) # Volume flow rate [m^3/s]
            Water_in.set_T(Tw_in) # Water entry temperature
    
            fan = [bool(int(sample[4])), bool(int(sample[5])), bool(int(sample[6])), bool(int(sample[7]))]
            HX.update_inputs(Water_in, Air_in, Water_out, Air_out, fan)
            print(HX.parameters)
            HX.solve()
            diff.append(Tw_out_exp-HX.Hxs[7].Tw_out)
    
    return diff

def optimize_Hx_ph2(x, C1, C2):
    
    
    C3 = x[0]
    C4 = x[1]


    with open("fake_sample.csv", "r") as fake_file:
        samples = csv.reader(fake_file, delimiter=";")
        diff = []
    
        for i, sample in enumerate(samples):
            if i==0:
                continue
            Tw_out_exp = float(sample[0])
            Tw_in = float(sample[1])
            Ta = float(sample[2])
            Qw = float(sample[3])
            fan = [bool(int(sample[4])), bool(int(sample[5])), bool(int(sample[6])), bool(int(sample[7]))]
            if fan != [True, True, True, True]:
                Hxs = Dry_cooler([C1, C2, C3, C4], Tw_in, Ta, Qw, fan)
                diff.append(Tw_out_exp-Hxs[7].Tw_out)
    
    return diff

def optimize_Hx_ph1(x):
    C1 = x[0]
    C2 = x[1]
    #C3 = x[2]
    #C4 = x[3]
    
    with open("fake_sample.csv", "r") as fake_file:
        samples = csv.reader(fake_file, delimiter=";")
        diff = []
    
        for i, sample in enumerate(samples):
            if i==0:
                continue
            Tw_out_exp = float(sample[0])
            Tw_in = float(sample[1])
            Ta = float(sample[2])
            Qw = float(sample[3])
            fan = [bool(int(sample[4])), bool(int(sample[5])), bool(int(sample[6])), bool(int(sample[7]))]
            if fan == [True, True, True, True]:
                Hxs = Dry_cooler([C1, C2, 1, 0], Tw_in, Ta, Qw, fan)
                diff.append(Tw_out_exp-Hxs[7].Tw_out)
    
    return diff

def fake_outputs(nb_sample, n1, n2):
    Qw_choice = np.arange(0.0008, 0.004, 0.0003)
    Tw_choice = np.arange(40, 65, 0.5)
    Ta_choice = np.arange(10, 20, 1)
    fan_choice = [0,1]
    
    
    with open("fake_sample.csv", "w") as fake_file:
        print("Tw_out ; Tw_in  ; Ta ; Qw ; fan1 ; fan2 ; fan3 ; fan4", file=fake_file)
        for i in range(nb_sample):
            Qw = random.choice(Qw_choice)
            Tw_in = random.choice(Tw_choice)
            Ta = random.choice(Ta_choice)
            fan = [1, 1, 1, random.choice(fan_choice)]
            Hxs = Dry_cooler([1, 0, 1, 0], Tw_in, Ta, Qw, fan)
            print("%.8f ; %.8f ; %.9f ; %.8f ; %s ; %s ; %s ; %s" %(Hxs[7].Tw_out + random.random()*n1+n2, Tw_in, Ta, Qw, fan[0], fan[1], fan[2], fan[3]), file=fake_file)
"Intput variables"
Qw = 0.0007 # Volume flow rate [m^3/s]
Tw = 65 # Water entry temperature
Ta = 20 # Air ambiant temperature
exch1_fan = True # If the fan is working or not
exch2_fan = True
exch3_fan = True
exch4_fan = False
#Hxs = exchanger_model([1,0,1,0], Tw, Ta, Qw, [True, True, True, False])
#fake_outputs(100, 3, -0.5)
C = optimize_HX(mode=1)
print(C)