# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 08:33:30 2024

@author: Elise
"""

from CoolProp.CoolProp import PropsSI

class Pump_cst_eff_is:
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
        #Citer tout les inputs
        #Supply state
        self.point_su = None
        #Output pressure
        self.point_ex = None
        self.N_pp = None
        
        # self.Required_inputs = [self.point_su.p, self.point_su.h, self.point_ex.p, self.N_pp]
        # self.check_calculable()

        "Outputs"
         #Citer tout les outputs
        
        "Parameter"
        self.epsilon_is = None # Isentropic efficiency
        self.V_s = None # Swept volume of the fluid
        self.epsilon_vol = None # Volumetric efficiency
        self.V = None # Volume of the pump

    def update_inputs(self, point_su, point_ex, N_pp):
        self.point_su = point_su
        self.point_ex = point_ex
        self.N_pp = N_pp
        
        # print(self.point_su.p, self.point_su.h, self.point_ex.p, self.N_pp)
        self.Required_inputs = [self.point_su.p, self.point_su.h, self.point_ex.p, self.N_pp]
        # print(self.Required_inputs)
        # print(self.Required_inputs)
        self.check_calculable()
        
    def check_calculable(self):
        # print("coucou")
        # print(self.Required_inputs)
        self.Required_inputs = [self.point_su.p, self.point_su.h, self.point_su.fluid, self.point_ex.p, self.N_pp]
        if all(Input is not None for Input in self.Required_inputs):   
            self.calculable = True
            
    def check_parametrized(self):
        if self.epsilon_is != None and self.V_s != None and self.epsilon_vol != None and self.V != None: 
            self.parametrized = True
            
        
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
                    print(f"Warning: Parameter '{key}' not found in the parameters.")

            self.check_parametrized()
               
        
    def solve(self):
        if self.calculable and self.parametrized:
            
            "Inputs"
            Fluid = self.point_su.fluid
            P_su = self.point_su.p
            h_su = self.point_su.h
            s_su = self.point_su.s
            rho_su = self.point_su.D
            P_ex = self.point_ex.p
            N_pp = self.N_pp

            h_min = PropsSI('H','P',5e4,'T',253.15,Fluid)
            h_max = PropsSI('H','P',4e6,'T',500,Fluid)

            #MODELLING PART
            if P_su < P_ex and N_pp > 0:      
            # If the external conditions are viable, we proceed to the modeling
                # print(P_ex, s_su, Fluid)
                try:
                    h_ex_s = PropsSI('H','P',P_ex,'S',s_su,Fluid)
                except:
                    h_ex_s = h_su
                m_dot = N_pp/60*self.V_s*self.epsilon_vol*rho_su
                W_dot = m_dot*(h_ex_s-h_su)/self.epsilon_is
                # Q_dot = AU*(T_su - T_amb)
                h_ex = h_su+W_dot/m_dot
                if h_ex > h_min and h_ex < h_max:
                    self.flag = 1
                else:
                    self.flag = -1
            
            self.point_ex.set_fluid(self.point_su.fluid)
            self.point_ex.set_m_dot(self.point_su.m_dot)
            self.point_ex.set_h(h_ex)
            # self.point_ex.set_p(P_ex)
            self.defined = True
            
        else:
            if self.calculable == False:
                print("Input of the component not completely known")
                
            if self.parametrized == False:
                print("Parameters of the component not completely known")
        

