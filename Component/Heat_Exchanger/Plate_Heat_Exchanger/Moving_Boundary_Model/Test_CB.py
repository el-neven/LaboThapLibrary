# from __future__ import division, print_function

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/LaboThapLibrary')

from Port.Mass_connector import Mass_connector
from HX_GeneralizedMovingBoundaries_Plate import Plate_HX_Geom_SWEP, Plate_HeatExchanger
from CoolProp.CoolProp import PropsSI  


"--------- 1) Evaporator ------------------------------------------------------------------------------------------"

"Refrigerant Su"
Evap_ref_su = Mass_connector()

# Set the fluid
Evap_ref_su.set_fluid('R1233zd(E)')

# Set other properties
Evap_ref_su.set_m_dot(0.4621)  # Example mass flow rate [kg/s] (from SWEP)
Evap_ref_su.set_T(59.98 + 273.15) # Example temperature [K] (from SWEP)
P_sat = PropsSI('P','T', 60+273.15,'Q', 0.5, Evap_ref_su.fluid) # Example Pressure [Pa]
Evap_ref_su.set_p(P_sat)  # Example Pressure [Pa]

"Water Su"
Evap_w_su = Mass_connector()
Evap_w_su.set_fluid("Water")

# Set other properties
Evap_w_su.set_m_dot(2.425)  # Example mass flow rate [kg/s]
Evap_w_su.set_T(70 + 273.15) # Example temperature [K]
Evap_w_su.set_p(2e5)  # Example Pressure [Pa]

"Pressure Drop"

DP_H_ON = False
DP_C_ON = False

"Heat Exchanger parameters"

n_disc = 15 # number of discretization

calc = 1 # flag to compute the HX
plot = 1 # flag to plot the HX fluids temperature profiles
print_flag = 1 # flag to print the HX results

"Geometry"

HX_Evap_geom = Plate_HX_Geom_SWEP()
HX_Evap_geom.set_parameters_SWEP("P200THx140/1P")

"Flow and htc correlation types parameters"

htc_type = "Correlation" # not user-defined but the relation are implemented in the code htc = heat transfer coefficient
flow_type = "Counter_flow"

"Heat Exchanger initiation and computation"

print("\n")
HX_Evap =  Plate_HeatExchanger()
HX_Evap.inputs(Evap_w_su, Evap_ref_su, wf_T='cold')
print(HX_Evap.calculable)
HX_Evap.set_parameters(**{
    'geom': HX_Evap_geom,
    'n_disc': n_disc
})
HX_Evap.solve()
print(HX_Evap.C_ex.h)
#HX_Evap.plot_objective_function()
HX_Evap.plot_cells()
# print(HX_Evap.Qmax)
# HX_Evap.plot_objective_function()

# "--------- 3) Condenser ------------------------------------------------------------------------------------------"

# "Cyclopentane Su"
# Cond_C5_su = Mass_connector()

# # Set the fluid
# Cond_C5_su.set_fluid("Cyclopentane")

# # Set other properties
# Cond_C5_su.set_m_dot(0.014)  # Example mass flow rate [kg/s]
# Cond_C5_su.set_T(139 + 273.15) # Example temperature [K]
# Cond_C5_su.set_p(0.9*1e5)  # Example Pressure [Pa]

# "Water Su"
# Cond_water_su = Mass_connector()

# # Set the fluid
# Cond_water_su.set_fluid("Water")

# # Set other properties
# Cond_water_su.set_m_dot(0.2)  # Example mass flow rate [kg/s]
# Cond_water_su.set_T(20 + 273.15) # Example temperature [K]
# Cond_water_su.set_p(5e5)  # Example Pressure [Pa]

# "Pressure Drop"

# DP_H_ON = True
# DP_C_ON = True

# "Heat Exchanger parameters"

# n_disc = 15 # number of discretization

# calc = 1 # flag to compute the HX
# plot = 0 # flag to plot the HX fluids temperature profiles
# print_flag = 1 # flag to print the HX results

# "Geometry"

# HX_Cond_geom = Plate_HX_Geom_SWEP()
# HX_Cond_geom.set_parameters_SWEP("B20Hx24/1P")

# "Flow and htc correlation types parameters"

# htc_type = "Correlation"
# flow_type = "Counter_flow"

# "Heat Exchanger"

# print("\n")
# HX_Cond, Cond_C5_ex, Cond_water_ex = Plate_HX("H",Cond_C5_su, Cond_water_su, HX_Cond_geom, htc_type, flow_type, n_disc, DP_H_ON, DP_C_ON, calc, plot, print_flag)

# M_HX_Cond = sum(HX_Cond.Mvec_h)

    