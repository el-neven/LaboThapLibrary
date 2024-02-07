# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 08:26:31 2024

@author: Elise
"""

from CoolProp.CoolProp import PropsSI

class Isenthalpic_valve:
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
        #Supply state
        self.su = None
        
        "Outputs"
        self.ex = None
        
        "Parameter"
        self.eta_is = None # Isentropic efficiency

        
    def check_calculable(self):
        if self.su != None and self.ex != None:
            if self.su.completely_known and self.ex.p != None:
                self.calculable = True
            
    def set_su(self, value):
        self.su = value
        self.check_calculable()
        
    def set_ex(self, value):
        self.ex = value
        self.check_calculable()
        
               
        
    def solve(self):
        if self.calculable:
            
            "Inputs"
            self.ex.set_fluid(self.su.fluid)
            self.ex.set_m_dot(self.su.m_dot)
            self.ex.set_h(self.su.h)
            self.defined = True
            
            
        else:
            if self.calculable == False:
                print("Input of the component not completely known")
                
            if self.parametrized == False:
                print("Parameters of the component not completely known")
        





