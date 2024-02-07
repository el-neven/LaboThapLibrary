# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 10:07:26 2024

@author: Elise
"""

import sys

sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/Python_Library')

import CoolProp
from CoolProp.CoolProp import PropsSI
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

from Component.Heat_Exchanger.Constant_Pinch_Point_Model.Simulation_Model import HX_Cst_PP
from Component.Valve.Isenthalpic_valve.Isenthalpic_valve import Isenthalpic_valve
from Component.Compressor.Constant_Isentopic_Efficiency.Simulation_Model import Compressor_cst_eff_is

from Port.Mass_connector import Mass_connector
#Import library
#from Saturation_curves import TS_curve_generator
#from Saturation_curves import PH_curve_generator


class Param_HP:
    def __init__(self, glide_cd, V_cp_swept,T_cd_htf_su, T_ev_htf_su, Fluid, RPM_cp, Pinch_ev, Pinch_cd, glide_ev, rv, rv_opt, eta_is_cp, eff_cp, Sc, Oh, dp_ev, dp_cd ):
        self.glide_cd = glide_cd #glide of temperature of condenser [K]
        self.V_cp_swept = V_cp_swept #swept volume of the compressor [m3]
        self.T_cd_htf_su = T_cd_htf_su #secondary fluid supply temperature of the condenser [°C]
        self.T_ev_htf_su = T_ev_htf_su #secondary fluid supply temperature of the evaporator [°C]
        self.Fluid = Fluid #working fluid
        self.RPM_cp = RPM_cp #compressor speed [RPM]
        self.Pinch_ev = Pinch_ev #pinch evaporator [K]
        self.Pinch_cd = Pinch_cd #pinch condenser [K]
        self.glide_ev = glide_ev #glide of temperature of the evaporator [K]
        self.rv = rv #volume ratio of the compressor
        self.rv_opt = rv_opt #rv_optim=0 means the volume ratio of the compressor is considered leading to under and overexpansion losses
        self.eta_is_cp = eta_is_cp #nominal maximum efficieny of the compressor
        self.eff_cp = eff_cp
        self.Sc = Sc #subcooling [K]
        self.Oh = Oh #Overheating [K]
        self.dp_ev = dp_ev #pressure drop in the evaporator line [bar]
        self.dp_cd = dp_cd #pressure drop in the condenser line [bar]
        

class HP:
    
    def __init__(self, Param):
        self.Param = Param

    def set_parameters(self, **kwargs):
        """
        Set parameters of the heat exchanger.

        Parameters
        ----------
        **kwargs : dict
            Key-value pairs representing parameters and their values.
            
            Example of call : heat_exchanger.set_parameters(D_o=0.05, t=0.002, H_core=2.0)
        """
        for key, value in kwargs.items():
            setattr(self, key, value)
    
    def Define(self): #Link components together
        #You need to put all of the connectors somewhere.

        "1. Condenseur"
        COND = HX_Cst_PP(Wf_su, Wf_ex, Sf_su, Sf_ex) #Quand même laisser les set??
        #Il faut savoir: tout les inputs,...

        "2. Isenthalpic valve"
        V = Isenthalpic_valve()

        "3. Evaporator"
        EVAP = HX_Cst_PP(Wf_su, Wf_ex, Sf_su, Sf_ex)

        "4. Compressor"
        COMP = Compressor_cst_eff_is()

    def System(self, x):
        
        Fluid = self.Param.Fluid
        #Iteration on the evaporation and condensation pressures
        self.P_ev, self.P_cd = x
        
        # Restrict the values of P_ev and P_cd within the desired range
        #Needs to find values to automize this!!!
        P_ev_max = PropsSI('P', 'T', self.Param.T_ev_htf_su, 'Q', 0, Fluid)
        #P_ev_min = ????
        P_cd_min = PropsSI('P', 'T', self.Param.T_cd_htf_su, 'Q', 0, Fluid)
        #P_cd_max = ?????
        
        self.P_ev = max(min(self.P_ev, P_ev_max), 10000)
        self.P_cd = max(min(self.P_cd, 1500000), P_cd_min-100000)
        print(self.P_ev, self.P_cd)
        #That doesn't move
        #--------------------------------------------------------------------------
        list_of_components = []

        "1. Condenser"
        #Définir tout les inputs du modèle ici??
        T_cd_sat = PropsSI('T', 'P', self.P_cd, 'Q', 0, Fluid)

        cd_wf_su = Mass_connector()
        cd_wf_ex = Mass_connector()
        cd_sf_su = Mass_connector()
        cd_sf_ex = Mass_connector()

        cd_wf_su.set_fluid(Fluid)
        cd_wf_su.set_m_dot(0.01) #Ici c'est inconnu donc normalement il doit passer ce composant ci en premier
        cd_wf_su.set_p(self.P_cd)
        cd_wf_su.set_T(T_cd_sat)

        cd_wf_ex.set_fluid(Fluid)
        cd_wf_ex.set_p(self.P_cd)
        cd_wf_ex.set_T(T_cd_sat-self.Sc)

        COND = HX_Cst_PP(cd_wf_su, cd_wf_ex, cd_sf_su, cd_sf_ex)

        #Supply
        #Condition at the compressor exit
        
        #Exhaust
        #T_cd_ex = T_cd_sat - self.Param.Sc
        #h_cd_ex = PropsSI('H', 'T', T_cd_ex, 'P', self.P_cd, Fluid)
        #s_cd_ex = PropsSI('S', 'T', T_cd_ex, 'P', self.P_cd, Fluid)
        #s_cd_1 = PropsSI('S', 'Q', 1, 'P', self.P_cd, Fluid) #For the graph
        #h_cd_1 = PropsSI('H', 'Q', 1, 'P', self.P_cd, Fluid) #For the graph
        #s_cd_0 = PropsSI('S', 'Q', 0, 'P', self.P_cd, Fluid) #For the graph
        #h_cd_0 = PropsSI('H', 'Q', 1, 'P', self.P_cd, Fluid) #For the graph
        #self.T_cd_ex = T_cd_ex
        #self.s_cd_ex = s_cd_ex
        #self.h_cd_ex = h_cd_ex
        #self.s_cd_1 = s_cd_1
        #self.h_cd_1 = h_cd_1
        #self.s_cd_0 = s_cd_0
        #self.h_cd_0 = h_cd_0
        
        "2. Isenthalpic valve"
        #valve su = cd ex
        #valve ex = evap su
        
        #--------------------------------------------------------------------------------
        "3. Evaporator"
        T_ev_sat = PropsSI('T', 'P', self.P_ev, 'Q', 0, Fluid)
        self.T_ev_sat = T_ev_sat
        
        #Supply
        h_ev_su = h_cd_ex; #Isenthalpic valve
        x_ev_su = PropsSI('Q', 'H', h_ev_su, 'P', self.P_ev, Fluid)
        s_ev_su = PropsSI('S', 'H', h_ev_su, 'P', self.P_ev, Fluid)
        self.s_ev_su = s_ev_su
        self.h_ev_su = h_ev_su
        
        #Exhaust
        T_ev_ex = T_ev_sat + self.Param.Oh
        rho_ev_ex = PropsSI('D', 'T', T_ev_ex,'P', self.P_ev, Fluid)
        mu_ev_ex = PropsSI('V', 'T', T_ev_ex,'P', self.P_ev, Fluid)
        h_ev_ex = PropsSI('H', 'T', T_ev_ex,'P', self.P_ev, Fluid)
        s_ev_ex = PropsSI('S', 'T', T_ev_ex,'P', self.P_ev, Fluid)
        s_ev_1 = PropsSI('S', 'Q', 1, 'P', self.P_ev, Fluid)
        h_ev_1 = PropsSI('H', 'Q', 1, 'P', self.P_ev, Fluid)
        self.s_ev_ex = s_ev_ex
        self.s_ev_1 = s_ev_1
        self.h_ev_ex = h_ev_ex
        self.h_ev_1 = h_ev_1
        self.T_ev_ex = T_ev_sat
        
        # Mass flow rate
        V_dot_cp = self.Param.RPM_cp/60*self.Param.V_cp_swept
        m_dot_wf = rho_ev_ex*V_dot_cp
        
        # Energy balance on the secondary fluid side
        Q_dot_ev = m_dot_wf*(h_ev_ex-h_ev_su)
        T_ev_htf_ex = self.Param.T_ev_htf_su-self.Param.glide_ev
        self.T_ev_htf_ex = T_ev_htf_ex
        cp_sf=4200 #Water (maybe put in parameters)
        m_dot_htf_ev = Q_dot_ev/(self.Param.glide_ev*cp_sf)

        #------------------------------------------------------------------------------
        "4. Compressor"
        
        #Supply
        P_cp_su = self.P_ev-self.Param.dp_ev
        self.P_cp_su = P_cp_su
        h_cp_su = PropsSI('H', 'T', T_ev_ex, 'P', P_cp_su, Fluid)
        s_cp_su = PropsSI('S', 'T', T_ev_ex, 'P', P_cp_su, Fluid)
        rho_cp_su = PropsSI('D','T', T_ev_ex, 'P', P_cp_su, Fluid)
        v_cp_su = 1/rho_cp_su
        Re_cp_su = rho_ev_ex*V_dot_cp/mu_ev_ex #Reynolds number
        
        #Exhaust
        P_cp_ex = self.P_cd+self.Param.dp_cd
        self.P_cp_ex = P_cp_ex
        h_cp_ex_is = PropsSI('H', 'S', s_cp_su, 'P', P_cp_ex, Fluid)
        
        if self.Param.rv_opt == 1:
            #If the rv is not taken into account
            h_cp_ex = h_cp_su + (h_cp_ex_is-h_cp_su)/self.Param.eff_cp
            self.h_cp_ex = h_cp_ex
        
        else:
            #If rv is taken into account -> compression in two parts:
            #Isentropic compression
            v_cp_in = v_cp_su/self.Param.rv
            h_cp_in = PropsSI('H', 'D', 1/v_cp_in, 'S', s_cp_su, Fluid)
            P_cp_in = PropsSI('P', 'D', 1/v_cp_in, 'S', s_cp_su, Fluid)
            w_in_v = h_cp_in-h_cp_su
            
            #Isochoric compression
            w_in_s = v_cp_in*(P_cp_ex-P_cp_in)
            
            #Total
            w_in = w_in_v+w_in_s
            h_cp_ex = h_cp_su + w_in/self.Param.eff_cp
            self.h_cp_ex = h_cp_ex
            
            W_dot_in_cp = m_dot_wf*w_in
        
        T_cp_ex = PropsSI('T', 'H', h_cp_ex, 'P', P_cp_ex, Fluid)
        s_cp_ex = PropsSI('S', 'H', h_cp_ex, 'P', P_cp_ex, Fluid)
        self.T_cp_ex = T_cp_ex
        self.s_cp_ex = s_cp_ex
        
        #Work made by the compressor (missing some losses to have the power at the shaft)
        self.W_dot_cp = m_dot_wf*(h_cp_ex-h_cp_su)
        
        #-------------------------------------------------------------------------------
        "Performances"
        self.COP = self.Q_dot_cd/self.W_dot_cp
        
        
        #--------------------------------------------------------------------------------
        "Pinch points"
        #Evaporator
        Pinch_ev = min(self.Param.T_ev_htf_su-T_ev_ex, T_ev_htf_ex-T_ev_sat)
        
        #Condensor
        #Revoir comment déterminer les pinch points
        h_cd_1 = PropsSI('H','Q',1,'P', self.P_cd, Fluid)
        self.Q_dot_cd = m_dot_wf*(h_cp_ex-h_cd_ex)
        Q_dot_cd_tpl = m_dot_wf*(h_cd_1-h_cd_ex)
        m_dot_htf_cd = self.Q_dot_cd/(cp_sf*self.Param.glide_cd)
        T_pinch_cd_1 = self.Param.T_cd_htf_su + Q_dot_cd_tpl/(m_dot_htf_cd*cp_sf)
        Pinch_cd = min(T_cd_sat-T_pinch_cd_1, T_cd_ex-self.Param.T_cd_htf_su)
        T_cd_htf_ex = self.Param.T_cd_htf_su + self.Param.glide_cd
        self.T_cd_htf_ex = T_cd_htf_ex
        
        "Residues"
        #Iteration until the pinch points calculated equals the ones set
        self.res_Pinch_ev = abs((self.Param.Pinch_ev-Pinch_ev)/Pinch_ev)
        self.res_Pinch_cd = abs((self.Param.Pinch_cd-Pinch_cd)/Pinch_cd)
        print(self.res_Pinch_ev, self.res_Pinch_cd)
        return self.res_Pinch_ev, self.res_Pinch_cd
        
    
    def solve_HP(self):
        #-------------------------------------------------------------------------
        P_ev_guess = PropsSI('P', 'T', self.Param.T_ev_htf_su, 'Q', 0, Fluid) - 300000
        P_cd_guess = PropsSI('P', 'T', self.Param.T_cd_htf_su, 'Q', 0, Fluid) + 430000
        #------------------------------------------------------------------------
        args = ()
        
        x = [P_ev_guess, P_cd_guess]
        
        result = fsolve(self.System, x, args=args)
    

#parameters found in on design model
glide_cd = 5 #glide of temperature of condenser [K]
glide_ev = 5 #glide of temperature of the evaporator [K]
V_cp_swept = 30/100**3 #swept volume of the compressor [m3]
T_cd_htf_su = 75+273.15 #secondary fluid supply temperature of the condenser [K]
T_ev_htf_su = 75+273.15 #secondary fluid supply temperature of the evaporator [K]
Fluid = 'R1233zd(E)' #working fluid
RPM_cp = 6000 #compressor speed [RPM]
Pinch_ev = 2 #pinch evaporator [K]
Pinch_cd = 2 #pinch condenser [K]
rv = 2.2 #volume ratio of the compressor
rv_opt = 0
eta_is_cp = 0.75 #nominal maximum isentropic efficieny of the compressor
eff_cp = 0.75 #Nominal maximum efficiency of the compressor
Sc = 5 #subcooling [K]
Oh = 5 #Over Heating [K]
dp_ev = 0.1*1**5 #pressure drop in the evaporator line [bar]
dp_cd = 0.1*1**5 #pressure drop in the condenser line [bar]

Param = Param_HP(glide_cd, V_cp_swept,T_cd_htf_su, T_ev_htf_su, Fluid, RPM_cp, Pinch_ev, Pinch_cd, glide_ev, rv, rv_opt, eta_is_cp, eff_cp, Sc, Oh, dp_ev, dp_cd )

Out = HP(Param)
Out.solve_HP()

#----------------------------------------------------------------------------------------
"Performances"

COP = Out.Q_dot_cd/Out.W_dot_cp

# #----------------------------------------------------------------------------------------
# "Plots"

# T_h = np.linspace(Out.T_cd_htf_ex, T_cd_htf_su, 100)
# s_h = np.linspace(Out.s_cp_ex, Out.s_cd_ex, 100)
# T_c = np.linspace(Out.T_ev_htf_ex, T_ev_htf_su, 100)
# s_c = np.linspace(Out.s_ev_su, Out.s_ev_ex, 100)

# #Create array with point of the cycle
# T_array = [Out.T_ev_sat, Out.T_ev_sat, Out.T_ev_ex, Out.T_cp_ex, Out.T_cd_sat, Out.T_cd_sat, Out.T_cd_ex, Out.T_ev_sat]
# s_array = [Out.s_ev_su, Out.s_ev_1, Out.s_ev_ex, Out.s_cp_ex,  Out.s_cd_1, Out.s_cd_0, Out.s_cd_ex, Out.s_ev_su]

# # Now, s_TS_curve and T_TS_curve contain the entropy-temperature data points
# TS = TS_curve_generator(Fluid)
# TS.TS_curve()
# s_TS_curve = TS.s_TS_curve
# T_TS_curve = TS.T_TS_curve

# # Plot the T-s diagram
# plt.figure()
# plt.plot(s_TS_curve, T_TS_curve, color='k', linestyle='--', label='T-s Curve')
# plt.plot(s_h, T_h, color='red', linestyle='-', label='Heat source')
# plt.plot(s_c, T_c, color='blue', linestyle='-', label='Heat sink')
# plt.xlabel('Entropy [J/kgK]')
# plt.ylabel('Temperature [K]')
# plt.title('Temperature-Entropy (T-s) Diagram')
# plt.grid(True)

# # Adding points from T_array and s_array
# plt.scatter(s_array, T_array, color='orange', marker='o', s=10, label='Cycle Points')

# # Join the points in T_array and s_array with lines
# plt.plot(s_array, T_array, color='orange', linestyle='-', linewidth=1)

# # Add a legend to distinguish the T-s curve, cycle points, ...
# plt.legend()

# # Enregistrez le graphique au format SVG
# plt.savefig('HP_mode_glide_25K_Ts.svg', format='svg')

# plt.show()

# #------------------------------------------------------------------------------------------------
# "ph"

# #Create array with point of the cycle
# P_array = [Out.P_ev, Out.P_ev, Out.P_cp_su, Out.P_cp_ex, Out.P_cd, Out.P_cd, Out.P_cd, Out.P_ev]
# h_array = [Out.h_ev_su, Out.h_ev_1, Out.h_ev_ex, Out.h_cp_ex,  Out.h_cd_1, Out.h_cd_0, Out.h_cd_ex, Out.h_ev_su]


# # Now, s_TS_curve and T_TS_curve contain the entropy-temperature data points
# PH = PH_curve_generator(Fluid)
# PH.PH_curve()
# P_PH_curve = PH.P_PH_curve
# h_PH_curve = PH.h_PH_curve

# # Plot the T-s diagram
# plt.figure()
# plt.plot(h_PH_curve, P_PH_curve, color='k', linestyle='--', label='P-h Curve')
# plt.xlabel('Enthalpy [J/kg]')
# plt.ylabel('Pressure [Pa]')
# plt.title('Pressure-Enthalpy (P-h) Diagram')
# plt.grid(True)

# # Adding points from T_array and s_array
# plt.scatter(h_array, P_array, color='orange', marker='o', s=10, label='Cycle Points')

# # Join the points in T_array and s_array with lines
# plt.plot(h_array, P_array, color='orange', linestyle='-', linewidth=1)

# # Add a legend to distinguish the T-s curve, cycle points, ...
# plt.legend()

# # Enregistrez le graphique au format SVG
# plt.savefig('HP_mode_glide_25K_Ph.svg', format='svg')

# plt.show()











































