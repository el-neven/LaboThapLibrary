# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 14:29:08 2024

@author: Elise
"""

from Port.Mass_connector import Mass_connector
from Isenthalpic_valve import Isenthalpic_valve

import numpy as np
import matplotlib.pyplot as plt

"Define class"
VR = Isenthalpic_valve()

"Set parameters"

"Define mass connector"
point_su = Mass_connector()
point_ex = Mass_connector()

"Inputs"
point_su.set_fluid('R1234yf')
point_su.set_m_dot(0.1) #[kg/s]
point_su.set_T(288.92) #[K]
point_su.set_p(400000) #[Pa]

point_ex.set_p(800000) #[Pa]


VR.set_su(point_su)
VR.set_ex(point_ex)

VR.solve()