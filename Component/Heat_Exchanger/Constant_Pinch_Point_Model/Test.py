# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 10:13:02 2024

@author: Elise
"""

from CoolProp.CoolProp import PropsSI

from Port.Mass_connector import Mass_connector
from Simulation_Model import HX_Cst_PP


#%%
"Condenseur"


"Mass connector"
# Working fluid
Wf_su = Mass_connector()
Wf_ex = Mass_connector()

# Secondary fluid
Sf_su = Mass_connector()
Sf_ex = Mass_connector()

"Inputs"
#Working fluids input
Wf_su.set_fluid('R1233zd(E)')
Wf_su.set_T(75+273.15)
Wf_su.set_p(500000)
Wf_su.set_m_dot(0.1)

T_cd_sat = PropsSI('T', 'P', 500000, 'Q', 0, 'R1233zd(E)')
Sc = 2
T_cd_ex = T_cd_sat - Sc
Wf_ex.set_T(T_cd_ex)

#Secondary fluids input
Sf_su.set_fluid('Water')
Sf_su.set_T(75+273.15)
Sf_su.set_m_dot(0.1)
print(Sf_su.m_dot)

"Definition of the model"



"Inputs of the model"
#Input can only come from ports
COND = HX_Cst_PP(Wf_su, Wf_ex, Sf_su, Sf_ex)
COND.set_parameters(glide=5, cp_sf=4600)

# "Mass connectors"
# # Working fluid
# Wf_su = Mass_connector()
# Wf_ex = Mass_connector()

# COND.set_su1(Wf_su)
# COND.set_ex1(Wf_ex)

# # Secondary fluid
# Sf_su = Mass_connector()
# Sf_ex = Mass_connector()

# COND.set_su2(Sf_su)
# COND.set_ex2(Sf_ex)

# "Inputs"
# #Working fluids input
# Wf_su.set_fluid('R1233zd(E)')
# Wf_su.set_T(75+273.15)
# Wf_su.set_p(500000)
# Wf_su.set_m_dot(0.1)

# #Secondary fluids input
# Sf_su.set_fluid('Water')
# Sf_su.set_T(75+273.15)

# "Parameters"
# COND.set_parameters(P_sat=500000, Sc=2, glide=5, cp_sf=4600)

COND.solve()

#%%
#Evaporator





