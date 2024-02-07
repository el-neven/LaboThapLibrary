

from CoolProp.CoolProp import PropsSI

from Model import HX_model_HP_CB

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/Python_Library')
from Port.Mass_connector import Mass_connector

# # EVAPORATOR

# #Test avec les données du POC1 en ORC mode
# #But: refaire le modèle du POC1

# #Mass connector
# point_su1 = Mass_connector()
# point_su2 = Mass_connector()

# point_ex1 = Mass_connector()
# point_ex2 = Mass_connector()

# P_in_wf = 4e5
# T_in_wf = PropsSI('T', 'P', P_in_wf, 'Q', 0, 'R1233zd(E)') -2 #Subcooled
# m_dot_sf = 0.0725032394688377

# T_out_wf = PropsSI('T', 'P', P_in_wf, 'Q', 1, 'R1233zd(E)') + 2 #Superheated

# T_in_sf = 72.5367743441399+273.15 #hot in the case of the evaporator
# m_dot_sf = 0.199615970833042

# # if PropsSI('H', 'T', T_in_wf, 'P', P_in_wf, 'R1233zd(E)') > PropsSI('H', 'P', P_in_wf, 'Q', 0, 'R1233zd(E)'):
# #     x_eva_in = 0
# # else:
# #     x_eva_in = -1

# point_su1.set_fluid('R1233zd(E)')
# point_su1.set_T(T_in_wf)
# # point_su1.set_p(P_in_wf)
# point_su1.set_m_dot(m_dot_sf)

# point_ex1.set_T(T_out_wf)

# point_su2.set_fluid('Water')
# point_su2.set_T(T_in_sf)
# point_su2.set_m_dot(m_dot_sf)
# point_su2.set_p(1e5)

# HX = HX_model_HP_CB()

# HX.set_parameters(**{
#     'HX_type': 'evaporator',
#     'HX_D': 0.02,
#     'HX_A': 2.38,
#     'min_pinch': 0.5,
#     'max_pinch': 20,
#     'glide_sf': 10,
#     'dT_sub_or_sup': 2
# })

# HX.update_inputs(point_su1, point_su2, point_ex1, point_ex2)

# HX.solve()

# print(point_ex1.h)

#CONDENSOR

#Test avec les données du POC1 en ORC mode
#But: refaire le modèle du POC1

#Mass connector
point_su1 = Mass_connector()
point_su2 = Mass_connector()

point_ex1 = Mass_connector()
point_ex2 = Mass_connector()

P_in_wf = 1.21923862687219e5
T_in_wf = PropsSI('T', 'P', P_in_wf, 'Q', 0, 'R1233zd(E)') + 5 # Superheated
m_dot_wf = 0.0725032394688377

T_out_wf = PropsSI('T', 'P', P_in_wf, 'Q', 1, 'R1233zd(E)') - 5 # Subcooled

T_in_sf = 14.3799641808547+273.15 #cold in the case of the condensor
m_dot_sf = 0.267194284887110

# if PropsSI('H', 'T', T_in_wf, 'P', P_in_wf, 'R1233zd(E)') > PropsSI('H', 'P', P_in_wf, 'Q', 0, 'R1233zd(E)'):
#     x_eva_in = 0
# else:
#     x_eva_in = -1

point_su1.set_fluid('R1233zd(E)')
point_su1.set_T(T_in_wf)
# point_su1.set_p(P_in_wf)
point_su1.set_m_dot(m_dot_wf)

point_ex1.set_T(T_out_wf)

point_su2.set_fluid('Water')
point_su2.set_T(T_in_sf)
point_su2.set_m_dot(m_dot_sf)
point_su2.set_p(1e5)

HX = HX_model_HP_CB()

HX.set_parameters(**{
    'HX_type': 'condenser',
    'HX_D': 0.02,
    'HX_A': 2.38,
    'min_pinch': 2,
    'max_pinch': 30,
    'glide_sf': 10,
    # 'dT_sub_or_sup': 2
})

HX.update_inputs(point_su1, point_su2, point_ex1, point_ex2)

HX.solve()

print(point_ex1.h)


# parity_plot(HP.p2,calc.p2); title('Evaporating pressure (bar)')
# parity_plot(HP.T2,calc.T2); title('Outlet temperature (°C)')
# parity_plot(HP.Q_eva./1000,calc.Q_eva); title('Thermal power (W)')
# parity_plot(HP.mc_sf,calc.mc_sf); title('Water flow rate (kg/s)')

# plot_figure(HP.m_wf,[eva.U_wf_sat,eva.U_wf_sub,eva.U_wf_sup],...
#     'WF mass flow rate (kg/s)','Heat transfer coefficient (W/K/m^2)','pointplot')
# legend('sat','sub','sup')

# plot_figure(HP.m_wf,[eva.A_sat,eva.A_sub,eva.A_sup],...
#     'WF mass flow rate (kg/s)','Heat transfer surface (m^2)','pointplot')
# legend('sat','sub','sup')
