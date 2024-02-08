"""
Supplemental code for paper:
I. Bell et al., "A Generalized Moving-Boundary Algorithm to Predict the Heat Transfer Rate of 
Counterflow Heat Exchangers for any Phase Configuration", Applied Thermal Engineering, 2014
"""

"""
Modification w/r to previous version:
    - Putting some order in the Objective Function "for" loops. Sparing some
    lines of code.
    - x_di_c correct calculation.
"""

# from __future__ import division, print_function

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/LaboThapLibrary')

from Port.Mass_connector import Mass_connector
from HX_GeneralizedMovingBoundaries_Plate import Plate_HX_Geom_SWEP, Plate_HeatExchanger
from CoolProp.CoolProp import PropsSI    

"--------- 1) Desuperheater ------------------------------------------------------------------------------------------"

"Cyclopentane Su"
DSH_C5_su = Mass_connector()

# Set the fluid
DSH_C5_su.set_fluid("Cyclopentane")

# Set other properties
DSH_C5_su.set_m_dot(0.014)  # Example mass flow rate [kg/s]
DSH_C5_su.set_T(205 + 273.15) # Example temperature [K]
DSH_C5_su.set_p(1*1e5)  # Example Pressure [Pa]

"Water Su"
DSH_water_su = Mass_connector()

# Set the fluid
DSH_water_su.set_fluid("Water")

# Set other properties
DSH_water_su.set_m_dot(0.08)  # Example mass flow rate [kg/s]
DSH_water_su.set_T(12 + 273.15) # Example temperature [K]
DSH_water_su.set_p(4e5)  # Example Pressure [Pa]

"Pressure Drop"

DP_H_ON = True
DP_C_ON = True

"Heat Exchanger parameters"

n_disc = 30 # number of discretization

calc = 1 # flag to compute the HX
plot = 0 # flag to plot the HX fluids temperature profiles
print_flag = 1 # flag to print the HX results

"Geometry"

HX_DSH_geom = Plate_HX_Geom_SWEP()
HX_DSH_geom.set_parameters_SWEP("B35TM0x10/1P")

"Flow and htc correlation types parameters"

htc_type = "Correlation"
flow_type = "Counter_flow"

"Heat Exchanger initiation and computation"

print("\n")

HX_DSH = Plate_HeatExchanger()

# HX_DSH, DSH_C5_ex, DSH_water_ex = Plate_HX("H",DSH_C5_su, DSH_water_su, HX_DSH_geom, htc_type, flow_type, n_disc, DP_H_ON, DP_C_ON, calc, plot, print_flag)

# M_HX_DSH = sum(HX_DSH.Mvec_h)

# "--------- 2) Evaporator ------------------------------------------------------------------------------------------"

# "Cyclopentane Su"
# Evap_C5_su = Mass_connector()

# # Set the fluid
# Evap_C5_su.set_fluid("Cyclopentane")

# # Set other properties
# Evap_C5_su.set_m_dot(0.014)  # Example mass flow rate [kg/s]
# Evap_C5_su.set_T(41 + 273.15) # Example temperature [K]
# Evap_C5_su.set_p(36*1e5)  # Example Pressure [Pa]

# "Oil Su"
# Evap_oil_su = Mass_connector()
# Evap_oil_su.set_fluid("INCOMP::T66")

# # Set other properties
# Evap_oil_su.set_m_dot(0.4)  # Example mass flow rate [kg/s]
# Evap_oil_su.set_T(243 + 273.15) # Example temperature [K]
# Evap_oil_su.set_p(5e5)  # Example Pressure [Pa]

# "Pressure Drop"

# DP_H_ON = True
# DP_C_ON = True

# "Heat Exchanger parameters"

# n_disc = 15 # number of discretization

# calc = 1 # flag to compute the HX
# plot = 0 # flag to plot the HX fluids temperature profiles
# print_flag = 1 # flag to print the HX results

# "Geometry"

# HX_Evap_geom = Plate_HX_Geom_SWEP()
# HX_Evap_geom.set_parameters_SWEP("B20Hx24/1P")

# "Flow and htc correlation types parameters"

# htc_type = "Correlation"
# flow_type = "Counter_flow"

# "Heat Exchanger initiation and computation"

# print("\n")
# HX_Evap, Evap_C5_ex, Evap_water_ex = Plate_HX("C",Evap_C5_su, Evap_oil_su, HX_Evap_geom, htc_type, flow_type, n_disc, DP_H_ON, DP_C_ON, calc, plot, print_flag)

# M_HX_Evap = sum(HX_Evap.Mvec_c)

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

    