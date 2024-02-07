# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 15:59:43 2024

@author: Elise
"""

from Component.Compressor.Constant_Isentopic_Efficiency.Simulation_Model import Compressor_cst_eff_is
from Component.Valve.Isenthalpic_valve.Isenthalpic_valve import Isenthalpic_valve
from Component.Heat_Exchanger.HX_Simple_model.HX_Simple_model import HX_simple

from Port.Mass_connector import Mass_connector
from Cycle import Cycle
from Cycle import Cycle_bis

"""-------------------
Components description
-------------------"""

list_of_components = []

""" Component 1 - COND """
COND = HX_simple()

COND.set_cp_sf(4.187)
COND.set_Delta_T_sf(25)

#list_of_components = list_of_components + [COND]


""" Component 2 - VR """
VR = Isenthalpic_valve()

# list_of_components = list_of_components + [VR]


""" Component 3 - EVAP """
EVAP = HX_simple()

EVAP.set_cp_sf(1.005)
EVAP.set_Delta_T_sf(-15)

list_of_components = list_of_components + [EVAP]
list_of_components = list_of_components + [VR]

""" Component 4 - COMP """
COMP = Compressor_cst_eff_is()

COMP.set_eta_is(0.8)

list_of_components = list_of_components + [COMP]

list_of_components = list_of_components + [COND]


"""--------------
Cycle description
--------------"""

cycle = Cycle(4, list_of_components)

"Inputs of the cycle"
#Condenser exhaust/ Valve supply
# cycle.point[0].set_fluid('R134a')
# cycle.point[0].set_T(45+273.15) #[K] -> Sinon over defined

#Valve exhaust/ Evaporator supply
# cycle.point[1].set_fluid('R134a')
cycle.point[1].set_p(800000) #[Pa]

#Comment prendre en compte les pertes de charges??? Là je les ai juste mise
#Evaporator exhaust/ Compressor supply
# cycle.point[2].set_fluid('R134a')
# cycle.point[2].set_p(800000-32000) #[Pa]

#Compressor exhaust/ Condenser supply
cycle.point[3].set_fluid('R134a')
cycle.point[3].set_m_dot(0.01)
cycle.point[3].set_T(105+273.15) #[Pa]
cycle.point[3].set_p(3000000) #[Pa]

"Secondary fluids"
water_cd_su = Mass_connector()
water_cd_ex = Mass_connector()

air_ev_su = Mass_connector()
air_ev_ex = Mass_connector()

water_cd_su.set_fluid('Water')
water_cd_su.set_m_dot(0.1)

water_cd_ex.set_fluid('Water')
#water_cd_ex.set_m_dot(0.1)

air_ev_su.set_fluid('Air')
air_ev_su.set_m_dot(0.1)

air_ev_ex.set_fluid('Air')
#air_ev_ex.set_m_dot(0.1)

"Component connection"
#Cycle mass connections
VR.set_su(cycle.point[0])
VR.set_ex(cycle.point[1])

EVAP.set_su1(cycle.point[1])
EVAP.set_ex1(cycle.point[2])

COMP.set_su(cycle.point[2])
COMP.set_ex(cycle.point[3])

COND.set_su1(cycle.point[3])
COND.set_ex1(cycle.point[0])

#Secondary fluid mass connections

EVAP.set_su2(air_ev_su)
EVAP.set_ex2(air_ev_ex)

COND.set_su2(water_cd_su)
COND.set_ex2(water_cd_su)

#Work and heat connections
cycle.solve()






















