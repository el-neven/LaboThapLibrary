# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 14:50:33 2024

@author: Elise
"""

from Port.Mass_connector import Mass_connector
from HX_Simple_model import HX_simple

import numpy as np
import matplotlib.pyplot as plt

# #-----------------------------------------------------------
# "Condenser example"
# COND = HX_simple()

# "Set parameters"

# COND.set_cp_sf(4.187)
# COND.set_Delta_T_sf(25)

# "Define mass connector"
# #Working Fluid
# point_su1 = Mass_connector()
# point_ex1 = Mass_connector()

# #Secondary fluid
# point_su2 = Mass_connector()
# point_ex2 = Mass_connector()

# "Inputs"
# point_su1.set_fluid('R410a')
# point_su1.set_m_dot(0.1) #[kg/s]
# point_su1.set_T(105+273.15) #[K]
# point_su1.set_p(3000000) #[Pa]

# point_su2.set_fluid('Water')
# point_su2.set_m_dot(0.1) #[Pa]

# COND.set_su1(point_su1)
# COND.set_su2(point_su2)

# COND.set_ex1(point_ex1)
# COND.set_ex2(point_ex2)

# COND.solve()

#-----------------------------------------------------------
"Evaporator example"
EVAP = HX_simple()

"Set parameters"

EVAP.set_cp_sf(1.005)
EVAP.set_Delta_T_sf(-15)

"Define mass connector"
#Working Fluid
point_su1 = Mass_connector()
point_ex1 = Mass_connector()

#Secondary fluid
point_su2 = Mass_connector()
point_ex2 = Mass_connector()

"Inputs"
point_su1.set_fluid('R410a')
point_su1.set_m_dot(0.1) #[kg/s]
point_su1.set_T(105+273.15) #[K]
point_su1.set_p(3000000) #[Pa]

point_su2.set_fluid('Air')
point_su2.set_m_dot(0.1) #[Pa]

EVAP.set_su1(point_su1)
EVAP.set_su2(point_su2)

EVAP.set_ex1(point_ex1)
EVAP.set_ex2(point_ex2)

EVAP.solve()
