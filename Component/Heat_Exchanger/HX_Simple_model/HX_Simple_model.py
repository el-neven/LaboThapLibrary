# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 14:32:39 2024

@author: Elise
"""


from CoolProp.CoolProp import PropsSI

class HX_simple:
    def __init__(self):
        """
        Parameters
        ----------
        /

        """

        self.calculable = False
        self.parametrized = False
        self.defined = False
        
        "Inputs"
        #Supply states
        self.su1 = None
        self.su2 = None
        
        "Outputs"
        self.ex1 = None
        self.ex2 = None

        
    def check_calculable(self):
        if self.su1 != None and self.su2 != None:
            if self.su1.completely_known and self.su2.m_dot != None and self.cp_sf != None and self.Delta_T_sf != None:
                self.calculable = True
            
    def set_su1(self, value):
        self.su1 = value
        self.check_calculable()
        
    def set_su2(self, value):
        self.su2 = value
        self.check_calculable()
        
    def set_ex1(self, value):
        self.ex1 = value
        self.check_calculable()
        
    def set_ex2(self, value):
        self.ex2 = value
        self.check_calculable()
        
    def set_cp_sf(self, value):
        self.cp_sf = value
        self.check_calculable()
        
    def set_Delta_T_sf(self, value):
        self.Delta_T_sf = value
        self.check_calculable()
        
               
        
    def solve(self):
        if self.calculable:
            
            "Inputs"
            h_su_wf = self.su1.h
            P_su_wf = self.su1.p
            T_su_wf = self.su1.T
            m_dot_wf = self.su1.m_dot

            m_dot_sf = self.su2.m_dot
            
            #Balance water side
            Q_dot = m_dot_sf*self.cp_sf*(self.Delta_T_sf)
            
            #Balance refrigerant side
            h_ex_wf = Q_dot/m_dot_wf + h_su_wf
            
            self.ex1.set_fluid(self.su1.fluid)
            self.ex1.set_m_dot(m_dot_wf)
            self.ex1.set_h(h_ex_wf)
            self.ex1.set_p(P_su_wf)
            self.defined = True
            
        else:
            if self.calculable == False:
                print("Input of the component not completely known")
                
            # if self.parametrized == False:
            #     print("Parameters of the component not completely known")






