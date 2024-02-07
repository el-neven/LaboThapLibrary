# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 16:10:24 2024

@author: Elise
"""

from CoolProp.CoolProp import PropsSI

class Compressor_cst_eff_is:
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
        #Output pressure
        self.ex = None
        
        "Outputs"
        
        "Parameter"
        self.eta_is = None # Isentropic efficiency

        
    def check_calculable(self):
        #Check si tout les inputs sont bien mit dedans!
        print(self.su.completely_known, self.ex)
        if self.su != None and self.ex != None:
            if self.su.completely_known and self.ex.p != None:
                self.calculable = True
        
    def check_parametrized(self):
        if self.eta_is != None:
            self.parametrized = True
            
    def set_su(self, value):
        self.su = value
        self.check_calculable()
        
    def set_ex(self, value):
        self.ex = value
        self.check_calculable()
        
    def set_eta_is(self, value):
        self.eta_is = value
        self.check_parametrized()
               
        
    def solve(self):
        if self.calculable and self.parametrized:
            
            "Inputs"
            P_ex = self.ex.p
            s_su = self.su.s
            h_su = self.su.h
            Fluid = self.su.fluid
            
            h_ex_is = PropsSI('H', 'P', P_ex, 'S', s_su, Fluid)
            h_ex = ((h_ex_is-h_su)/self.eta_is) + h_su
            
            self.ex.set_fluid(self.su.fluid)
            self.ex.set_m_dot(self.su.m_dot)
            self.ex.set_h(h_ex)
            self.defined = True
            
            
        else:
            if self.calculable == False:
                print("Input of the component not completely known")
                
            if self.parametrized == False:
                print("Parameters of the component not completely known")
        





