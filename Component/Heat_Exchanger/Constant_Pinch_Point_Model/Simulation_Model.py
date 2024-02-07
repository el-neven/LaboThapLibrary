# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 10:12:54 2024

@author: Elise
"""

from CoolProp.CoolProp import PropsSI

#Constant Pinch point Model of a Heat Exchanger

class HX_Cst_PP:
    def __init__(self, WF_su, WF_ex, SF_su, SF_ex):
        """
        Initialize the Simulation_Model with working fluid and secondary fluid at start-up and exhaust.
        """
        self.defined = False
        self.calculable = False

        # Working fluid properties
        self.WF_su = WF_su
        self.WF_ex = WF_ex
        self.T_su_wf = WF_su.T
        self.h_su_wf = WF_su.h
        self.P_sat = WF_su.p
        self.m_dot_wf = WF_su.m_dot
        self.T_ex_wf = WF_ex.T

        # Secondary fluid properties
        self.SF_su = SF_su
        self.SF_ex = SF_ex
        self.T_su_sf = SF_su.T
        self.m_dot_sf = SF_su.m_dot

        # Required inputs for the model
        self.Required_inputs = [self.T_su_wf, self.P_sat, self.m_dot_wf, self.T_ex_wf, self.T_su_sf, self.m_dot_sf, self.T_ex_wf]
        
        # Additional parameters
        self.glide = None  # Add glide as an instance variable
        self.cp_sf = None  # Add cp_sf as an instance variable
        
        self.check_calculable()
        
    def check_calculable(self):
    # if self.su1 != None and self.su2 != None:
        if all(Input is not None for Input in self.Required_inputs):
            self.calculable = True
        
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
                if hasattr(self, key):
                    setattr(self, key, value)
                else:
                    print(f"Warning: Parameter '{key}' not found in the heat exchanger.")
        
        
    def solve(self):
        if self.calculable:
            
            "Supply state"
            T_su_wf = self.T_su_wf
            h_su_wf = self.h_su_wf
            m_dot_wf = self.m_dot_wf
            WF = self.WF_su.fluid
            
            T_su_sf = self.T_su_sf
            m_dot_sf = self.m_dot_sf
            SF = self.SF_su.fluid
            
            self.T_sat_wf = PropsSI('T', 'P', self.P_sat, 'Q', 0, WF)
            
            "Echaust state"
            h_ex_wf = PropsSI('H', 'T', self.T_ex_wf, 'P', self.P_sat, WF)
            h_1 = 1900 #○PropsSI('H', 'Q', 1,' P', self.P_sat, WF)
            
            "Mass flow rate calculation"
            # print(h_su_wf, h_ex_wf)
            Q_dot = m_dot_wf*(h_ex_wf-h_su_wf)
            m_dot_sf = Q_dot/(self.cp_sf*self.glide)
            
            "Pinch point calculation"
            Q_dot_cd_tpl = m_dot_wf*(h_1-h_ex_wf)
            print(m_dot_sf, self.cp_sf)
            T_pinch_cd_1 = T_su_wf + Q_dot_cd_tpl/(m_dot_sf*self.cp_sf)
            print(m_dot_sf, self.cp_sf)
            self.Pinch_cd = min(self.T_sat_wf-T_pinch_cd_1, self.T_ex_wf-T_su_sf)
            
            T_ex_sf = T_su_sf + self.glide
            #Balance water side
            # Q_dot = m_dot_sf*self.cp_sf*(self.Delta_T_sf)
            
            self.WF_ex.set_fluid(self.WF_su)
            self.WF_ex.set_m_dot(m_dot_wf)
            self.WF_ex.set_p(self.P_sat)
            #self.WF_ex.set_T(self.T_ex_wf)
            self.defined = True
            
        else:
            if self.calculable == False:
                print("Input of the component not completely known")
                
            # if self.parametrized == False:
            #     print("Parameters of the component not completely known")






