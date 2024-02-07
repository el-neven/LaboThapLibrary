# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 08:33:30 2024

@author: Elise
"""

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/Python_Library')

from Port.Mass_connector import Mass_connector
from Component.Pump.Constant_Efficiency.Model import Pump_cst_eff_is

import numpy as np
import matplotlib.pyplot as plt

"Define mass connectors"
point_su = Mass_connector()
point_ex = Mass_connector()

"Define work connector"
#Pas encore fait (Attention vraiment important de faire!!, Ici on va juste le mettre dans les inputs)
N_pp = 1500


"Inputs"
point_su.set_fluid('R245fa')
# point_su.set_m_dot(0.1) #[kg/s] #Pas besoin de définir le débit massique, il est calculé par le modèle
point_su.set_h(2.6676e+05) #[K] 
point_su.set_p(4.0001e5) #[Pa]

point_ex.set_p(3.6510e+06*0.99) #[Pa]

"Define class"
PUMP = Pump_cst_eff_is(point_su, point_ex, N_pp)

"Set parameters"
PUMP.set_parameters(epsilon_is = 0.5, epsilon_vol = 0.8, V_s = 1e-6, V = 1.4e-3)


PUMP.solve()