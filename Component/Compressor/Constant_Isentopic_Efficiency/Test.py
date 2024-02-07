# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 08:33:30 2024

@author: Elise
"""

from Port.Mass_connector import Mass_connector
from Simulation_Model import Compressor_cst_eff_is

import numpy as np
import matplotlib.pyplot as plt

"Define class"
COMP = Compressor_cst_eff_is()

"Set parameters"
COMP.set_eta_is(0.8)

"Define mass connector"
point_su = Mass_connector()
point_ex = Mass_connector()

"Inputs"
point_su.set_fluid('R1234yf')
point_su.set_m_dot(0.1) #[kg/s]
point_su.set_T(288.92) #[K]
point_su.set_p(400000) #[Pa]

point_ex.set_p(800000) #[Pa]


COMP.set_su(point_su)
COMP.set_ex(point_ex)

COMP.solve()