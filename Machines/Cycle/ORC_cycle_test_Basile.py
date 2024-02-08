# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 11:47:52 2024
    
@author: bchaudoir
"""


# import CoolProp
from CoolProp.CoolProp import PropsSI
# from CoolProp.Plots import PropertyPlot
# import warnings

import matplotlib.pyplot as plt
import numpy as np

# import matplotlib.pyplot as plt
# import numpy as np

# from scipy.optimize import fsolve

# from Modules.EffExpander_V3 import EffExpander
# from Modules.ThermoState_221211 import ThermoState as ThermoState
# from Modules.ThermoState_Mix import ThermoState_Mix
# from Modules.HX_LMTD_calcUA_210610 import HX_LMTD_calcUA

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/LaboThapLibrary')

from Port.Mass_connector import Mass_connector

from Component.Compressor.Constant_Isentopic_Efficiency.Simulation_Model import Compressor_cst_eff_is
from Component.Valve.Isenthalpic_valve.Isenthalpic_valve import Isenthalpic_valve
from Component.Heat_Exchanger.HX_Simple_model.HX_Simple_model import HX_simple
from Component.Heat_Exchanger.Plate_Heat_Exchanger.Moving_Boundary_Model.HX_GeneralizedMovingBoundaries_Plate import Plate_HX_Geom_SWEP, Plate_HeatExchanger

from Component.Expander.Semi_Empirical_Model_VL.Simulation_Model import Expander_SE
from Component.Heat_Exchanger.Simple_Moving_Boundaries_Model.Model import HX_model_HP_CB
from Component.Pump.Constant_Efficiency.Model import Pump_cst_eff_is

from Port.Mass_connector import Mass_connector
from Port.Heat_connector import Heat_connector

from Cycle import Cycle

from scipy.optimize import fsolve
from scipy.optimize import root
from scipy.optimize import minimize
from scipy.optimize import least_squares

# from Modules.oil_point import Oil_Point_on_cycle
# from HX_GeneralizedMovingBoundaries_Plate import Plate_HX_Geom_SWEP, Plate_HX

#%%

class Thermal_System:
    #-------------------------------------------------------------------------
    class Cycles():
        pass
    
    class Components():
        pass
    
    class Components_chain_element():
        def __init__(self, name, nb_input, nb_output, cycle_names):
            self.name = name
            self.previous = [None]*nb_input # will be a thermodynamical point
            self.next = [None]*nb_output # will be a thermodynamical point
            self.cycles = cycle_names # Cycle(s) that the component is a part of (s in case of an HTX) 
            
    class Sources():
        pass
    
    class Source_chain_element():
        def __init__(self, name):
            self.name = name
            self.next = [None] # will be a thermodynamical point
            
    class Sinks():
        pass
    
    class Sink_chain_element():
        def __init__(self, name):
            self.name = name
            self.previous = [None] # will be a thermodynamical point
            
    class Params():
        pass
    #-------------------------------------------------------------------------
    
    def __init__(self, fluid, h_fluid, c_fluid):
        """
        Inputs:
        -------
        fluid : ORC working fluid
        h_fluid : Hot source fluid/mixture
        c_fluid : Cold source fluid/mixture
        --------------------------------------------------
        
        An upgrade would be to enter as input an array of fluids in order to have something else than 3 subcycles in the system
        """
        
        self.parametrized = False
        
        self.finished_it = False
        self.start_key = None
        
        self.Params.fluid = fluid # Working Fluid
        self.Params.h_fluid = h_fluid # Hot source fluid
        self.Params.c_fluid = c_fluid # Cold source fluid
        
        self.Components = self.Components()
        self.Components.dict = {} # Creates a dictionnary that has as keys the name of the components and as values their models
        self.Components.chains = {fluid : {}, h_fluid : {}, c_fluid : {}} # Create chains that will contain for each cycles the components and thermodynamical points -> Will be iterated on
        
        self.Sources = self.Sources()
        self.Sources.dict = {} # Creates a dictionnary to contain all sources (no model or inputs will be passed to them, sources and sinks just ensure that something is there without having a loop)

        self.Sinks = self.Sinks()
        self.Sinks.dict = {} # Creates a dictionnary to contain all sources

        self.Cycles = {fluid : {}, h_fluid : {}, c_fluid : {}} # Will contain only the thermodynamical points related to the cycles 

    #-------------------------------------------------------------------------
#%%

    def add_component(self, component, name, nb_input, nb_output, cycle_names):
        if name in self.Components.dict:
            print(f"Component '{name}' already added.")
            return
        
        # Add the component to the dictionnary
        self.Components.dict[name] = component
        
        # Create the component chain element
        new_chain_component = self.Components_chain_element(name, nb_input, nb_output, cycle_names)
        
        for cycle_name in cycle_names:
            self.Components.chains[cycle_name][name] = new_chain_component
            
    def link_components(self, name1, output_1, name2, input_2, cycle_name):
                
        if self.Components.chains[cycle_name][name1].next[output_1-1] != None or self.Components.chains[cycle_name][name2].previous[input_2-1] != None:
            print('One of ',name1, 'or', name2, 'is already linked at the specified port' )
            return
        else:
            # Create a thermodynamical point and list it in the cycles and component.chains sub_objects
            new_point_key = cycle_name[0] + str(len(self.Cycles[cycle_name]) + 1)
            self.Cycles[cycle_name][new_point_key] = Mass_connector()
            self.Components.chains[cycle_name][new_point_key] = self.Cycles[cycle_name][new_point_key] 
            
            # Link the component 1 to the thermodynamical point
            self.Components.chains[cycle_name][name1].next[output_1-1] = self.Cycles[cycle_name][new_point_key]
            self.Cycles[cycle_name][new_point_key].previous = self.Components.chains[cycle_name][name1]
            self.Components.dict[name1].point_ex[output_1-1] = self.Cycles[cycle_name][new_point_key]
            
            # Link the thermodynamical point to the component 2
            self.Components.chains[cycle_name][name2].previous[input_2-1] = self.Cycles[cycle_name][new_point_key]
            self.Cycles[cycle_name][new_point_key].next = self.Components.chains[cycle_name][name2]
            self.Components.dict[name2].point_su[input_2-1] = self.Cycles[cycle_name][new_point_key]
            
            # Set the fluid of the connector
            self.Cycles[cycle_name][new_point_key].set_fluid(cycle_name)        
    
                    # oil_source,"Oil_source_1",h_fluid,'Evaporator',2
    def add_source(self, source, name, cycle_name, name_next_comp, input_port):
        if name in self.Sources.dict:
            print(f"Source '{name}' already added.")
            return
        
        # Put the source in the source dictionnary
        self.Sources.dict[name] = source
        component = self.Components.chains[cycle_name][name_next_comp]

        # Create the source chain element
        Source_chain = self.Source_chain_element(name)
        
        # Create a thermodynamical point and add it to the subcycle
        new_point_key = cycle_name[0] + str(len(self.Cycles[cycle_name]) + 1)                    
        self.Cycles[cycle_name][new_point_key] = Mass_connector()
        
        # Link the source to the thermodynamical point
        self.Components.chains[cycle_name][new_point_key] = self.Cycles[cycle_name][new_point_key] 
        Source_chain.next = self.Cycles[cycle_name][new_point_key]
        self.Cycles[cycle_name][new_point_key].previous = Source_chain
        
        # Link the thermodynamical point to the component after the source
        component.previous[input_port-1] = self.Cycles[cycle_name][new_point_key]
        self.Cycles[cycle_name][new_point_key].next = component
        self.Components.dict[name_next_comp].point_su[input_port-1] = self.Cycles[cycle_name][new_point_key]

        # Set the fluid of the connector
        self.Cycles[cycle_name][new_point_key].set_fluid(cycle_name)
        
    def add_sink(self, sink, name, cycle_name, name_previous_comp,output_port):
        if name in self.Sinks.dict:
            print(f"Source '{name}' already added.")
            return
        
        # Put the sink in the source dictionnary
        self.Sinks.dict[name] = sink
        component = self.Components.chains[cycle_name][name_previous_comp]
        
        # Create the source chain element   
        Sink_chain = self.Sink_chain_element(name)
        
        # Create a thermodynamical point and add it to the subcycle
        new_point_key = cycle_name[0] + str(len(self.Cycles[cycle_name]) + 1)    
        self.Cycles[cycle_name][new_point_key] = Mass_connector()
        
        # Link the thermodynamical point to the sink
        self.Components.chains[cycle_name][new_point_key] = self.Cycles[cycle_name][new_point_key] 
        Sink_chain.previous = self.Cycles[cycle_name][new_point_key]
        self.Cycles[cycle_name][new_point_key].next = Sink_chain
        
        # Link the the component before the sink to the thermodynamical point
        component.next[output_port-1] = self.Cycles[cycle_name][new_point_key]
        self.Cycles[cycle_name][new_point_key].previous = component
        self.Components.dict[name_previous_comp].point_ex[output_port-1] = self.Cycles[cycle_name][new_point_key]
        
        # Set the fluid of the connector
        self.Cycles[cycle_name][new_point_key].set_fluid(cycle_name)        

#%%

    def check_parametrized(self):
        """
        Check for all components models if all the required parameters have been passed
        """        
        flag_not_parametrized = 0 # =1 if at least one component is not fully parametrized
        
        for comp in self.Components.dict:
            if self.Components.dict[comp].parametrized:
                pass
            else:
                self.defined = False
                flag_not_parametrized = 1
                print(f"Error: Parameters of {comp} not completely known")
        
        if flag_not_parametrized == 0:
            self.parametrized = True
            print("The cycle components are fully parametrized")

    def check_connections(self):
        """
        Just checks if all components have something in their inputs / outputs
        Possible upgrade : also check that all components are part of a cycle or are between a source and a sink

        """
        
        for key_chain in self.Components.chains:
            iter_chain = self.Components.chains[key_chain]
            
            for key_iter_component in iter_chain:
                iter_component = iter_chain[key_iter_component]
            
                if (not (isinstance(iter_component, self.Sink_chain_element) or isinstance(iter_component, Mass_connector))):
                    
                    for next_comp in iter_component.next:
                        if next_comp == None:
                            print("There is a missing 'next' connection for the component : ", iter_component.name)
                            self.Params.is_cycle = False
                            return False
                
                if (not (isinstance(iter_component, self.Source_chain_element) or isinstance(iter_component, Mass_connector))):
                    for prev_comp in iter_component.previous:
                        if prev_comp == None:
                            print("There is a missing 'previous' connection for the component : ", iter_component.name)
                            self.Params.is_cycle = False
                            return False
        
        self.Params.is_cycle = True
        return True

#%% 

    def set_su(self, chain_name, comp_name, port, variable, value):
        
        component_chain_elem = self.Components.chains[chain_name][comp_name]
        connector = component_chain_elem.previous[port]
        
        if variable == 'T':
            connector.set_T(value)            
        elif variable == 'P':
            connector.set_p(value)
        elif variable == 'M_dot':
            for connector in self.Cycles[chain_name]:
                self.Cycles[chain_name][connector].set_m_dot(value)
            
        return
    
    def set_ex(self, chain_name, comp_name, port, variable, value):
        
        component_chain_elem = self.Components.chains[chain_name][comp_name]
        connector = component_chain_elem.next[port]
        
        if variable == 'T':
            connector.set_T(value)            
        elif variable == 'P':
            connector.set_p(value)
        elif variable == 'M_dot':
            for key in self.Cycles[chain_name]:
                connector = self.Cycles[chain_name][key]
                connector.set_m_dot(value)
            
        return

#%% 
    def solve_chain_element(self, key_chain_elem, cycle, start_key):
        print("START KEY : ", start_key)
        
        if self.finished_it == False: # Check if all components were already computed for this iteration
            print(key_chain_elem)
            chain_element = self.Components.chains[cycle][key_chain_elem]
            
            if isinstance(chain_element, Mass_connector): # We are at a Mass connector
                print("Mass connector")
                
                if isinstance(chain_element.next, self.Sink_chain_element):
                    print("Next is a sink")
                    print("----------------")
                    return
                
                if chain_element.state_known: # State is known, try to solve the next component
                    print("State known !!!!")
                    
                    if self.start_key == None: # Real start of the cycle iteration -> saved
                        self.start_key = key_chain_elem
                        print("Set of the start key")

                    print("----------------")

                    next_comp = chain_element.next # Shall be a component
                    
                    # Right now I am using the keys to know beter what is happening, it might be possible to only pass the chain_element instead of the key to the method (that would be more efficient)
                    # Iterate through the items and find the key of the next component
                    for key, value in self.Components.chains[cycle].items():
                        if value == next_comp:
                            next_key = key
                            break
                        
                    if self.start_key == None:
                        self.solve_chain_element(next_key,cycle,start_key)
                    else:
                        self.solve_chain_element(next_key,cycle,self.start_key)
                    
                else: # State is not known, try to look directly for the next mass connector
                    print("State not known")
                    
                    next_key = None
                    for next_point in chain_element.next.next:
                        if next_key != None:
                            break
                        # Iterate through the items and find the key
                        for key, value in self.Components.chains[cycle].items():
                            if value == next_point:
                                next_key = key
                                break
                            
                    print("----------------")
                    
                    if self.start_key == None:
                        if next_key != start_key:
                            self.solve_chain_element(next_key,cycle,start_key)
                        else:
                            print("Back to start key")
                            print("!!!!!!!!!!!!!!!!!")
                            return
                    else:
                        if next_key != self.start_key:
                            self.solve_chain_element(next_key,cycle,self.start_key)   
                        else:
                            print("Back to start key")
                            print("!!!!!!!!!!!!!!!!!")
                            return
                        
            else: # We are at a Component
            
                next_key = None
                for next_point in chain_element.next:
                    if next_key != None:
                        break
                    # Iterate through the items and find the key
                    for key, value in self.Components.chains[cycle].items():
                        if value == next_point:
                            next_key = key
                            break  
                
                if next_key != self.start_key:
                    comp_model = self.Components.dict[key_chain_elem]
                    
                    if comp_model.defined == False:
                        comp_model.solve()
                        print("Component Solved")
                        print("????????????????")
                        
                        for i in range(len(comp_model.point_ex)):
                            
                            chain_element.next[i].set_p(comp_model.point_ex[i].p)                        
                            chain_element.next[i].set_h(comp_model.point_ex[i].h)  
                            
                            # print("- - - - - - - - - - ")
        
                            # chain_element.next[i].print_resume()
                            
                            # print("- - - - - - - - - - ")
                        
                        self.solve_chain_element(next_key,cycle,start_key)
                    else: 
                        print("Component already known")
                        print("----------")
                        
                else: # Dernier point à calculer
                    comp_model = self.Components.dict[key_chain_elem]
                    comp_model.solve()

                    for i in range(len(comp_model.point_ex)):
                        
                        chain_element.next[i].set_p(comp_model.point_ex[i].p)                        
                        chain_element.next[i].set_h(comp_model.point_ex[i].h)  
                        
                    print("Back to start key")
                    print("----------")
                    
                    self.finished_it = True
        
        else:
            print("All components are computed")
        
        return

    def solve(self):
        if self.parametrized == True:
            res = 1e10
            n_it_max = 1
            
            tol = 1
            n_it = 0
            while res > tol and n_it < n_it_max:
                print("--------------------------------------\nIteration :",n_it)
                for cycle in self.Cycles:
                    start_key = next(iter(self.Cycles[cycle]))
                    self.solve_chain_element(start_key, cycle, start_key)
                    
                    #Search the next component calculable
                    # if not self.compo_list[i].calculable:
                    #     self.compo_list[i].check_calculable()

                    # if self.compo_list[i].calculable and not self.compo_list[i].defined: #If the component is calculable and not defined then we solve it
                    #     self.compo_list[i].solve()
                    #     n_solved +=1

                    # else:
                    #     pass
                n_it +=1
            if res < tol and n_it < n_it_max:
                print("Cycle solved")
                
            else:
                print("Error: Cycle couldn't be solved")

        else:
            print("Error: Parameters of all components not completely known")

#%%
    
    def plot_cycle(self, fluid_name):
        fig, ax = plt.subplots()
        ax.axis('off')  # Turn off the axis
    
        components = self.Components.chains[fluid_name]
        coordinates = {}
    
        for i, (comp_name, component) in enumerate(components.items()):
            x = i * 100
            y = 0
            coordinates[comp_name] = (x, y)
            if isinstance(component, Mass_connector):
                ax.plot(x, y, 'o', markersize=12, color='black', fillstyle='none', markeredgecolor='white')  # Represent Mass_connector with circles (black border)
            else:
                ax.text(x, y, comp_name, ha='center', va='center', fontsize=8, bbox=dict(facecolor='white', edgecolor='black'))
    
        for i, (point_name, point) in enumerate(self.Cycles[fluid_name].items()):
            x, y = coordinates[point_name]
            ax.plot(x, y, 'o', markersize=12, color='white', markeredgecolor='black')  # Use markeredgecolor to represent the white border
            ax.text(x, y, point_name, ha='center', va='center', fontsize=8)
    
        plt.show()
        
    def show_cycle_info(self,wf):
        
        print("-----------------------")
        print("-----------------------")
        print("Cycle Results")
        print("-----------------------")
                
        for key in self.Cycles[wf]:
            print(key)
            
            point = self.Cycles[wf][key]
            
            print("T :", point.T)
            print("P :", point.p)
            print("H :", point.h)
            print("---------")
            
            plt.plot(point.h, point.p, marker='o', label=key, color = 'green')
        
        W_dot_pump = (self.Components.chains[wf]["Pump"].next[0].h - self.Components.chains[wf]["Pump"].previous[0].h)*self.Components.chains[wf]["Pump"].previous[0].m_dot
        W_dot_exp = (self.Components.chains[wf]["Expander"].previous[0].h - self.Components.chains[wf]["Expander"].next[0].h)*self.Components.chains[wf]["Expander"].previous[0].m_dot
        Q_dot_evap = (self.Components.chains[wf]["Evaporator"].next[0].h - self.Components.chains[wf]["Evaporator"].previous[0].h)*self.Components.chains[wf]["Evaporator"].previous[0].m_dot
        Q_dot_cond = (self.Components.chains[wf]["Condenser"].previous[0].h - self.Components.chains[wf]["Condenser"].next[0].h)*self.Components.chains[wf]["Condenser"].previous[0].m_dot
        
        print("W_dot_pump :",W_dot_pump," [W]")
        print("W_dot_exp :",W_dot_exp," [W]")
        print("Q_dot_evap :",Q_dot_evap," [W]")
        print("Q_dot_cond :",Q_dot_cond," [W]")

        print("-----------------------")

        # Saturation curve
        P_sat = np.linspace(1e5,45e5,1001)
        
        h_sat1 = np.zeros(len(P_sat))
        h_sat2 = np.zeros(len(P_sat))
        
        for i in range(len(h_sat1)):
            h_sat1[i] = PropsSI('H','P',P_sat[i],'Q',0,'Cyclopentane')
        
        for i in range(len(h_sat2)):
            h_sat2[i] = PropsSI('H','P',P_sat[i],'Q',1,'Cyclopentane')
        
        plt.plot(h_sat1, P_sat, label=key, color = 'black')
        plt.plot(h_sat2, P_sat, label=key, color = 'black')

        # Add labels and legend
        plt.xlabel('Enthalpy (h)')
        plt.ylabel('Pressure (p)')
        plt.title('P-h Diagram')
        plt.grid()
        plt.show()
        
        return
    
    # def circular_distribution(self, num_points, radius):
    #     angles = np.linspace(0, 2*np.pi, num_points, endpoint=False)
    #     x = radius * np.cos(angles)
    #     y = radius * np.sin(angles)
    #     return x, y
    
    # def plot_cycle(self, fluid_name):
    #     fig, ax = plt.subplots()
    #     ax.axis('off')  # Turn off the axis
        
    #     components = self.Components.chains[fluid_name]
    #     print(components.keys())
        
    #     num_components = len(components)
        
    #     x, y = self.circular_distribution(num_components, radius=100)  # Adjust the radius as needed
                
    #     coordinates = {comp_name: (x[i], y[i]) for i, comp_name in enumerate(components)}
        
    #     print(coordinates)
        
    #     # Inside the loop where you're plotting lines between components
    #     for comp_name, component in components.items():
    #         if isinstance(component, Mass_connector):
    #             continue  # Skip drawing lines for Mass_connector
            
    #         x1, y1 = coordinates[comp_name]
    #         for next_point in component.next:
    #             if isinstance(next_point, Mass_connector):
    #                 next_point_key = list(next_point.previous.keys())[0]  # Get the key for the Mass_connector
    #                 x2, y2 = coordinates[next_point_key]
    #                 ax.plot([x1, x2], [y1, y2], 'k-', lw=1)  # Draw lines between components

        
    #     for comp_name, (x, y) in coordinates.items():
    #         if isinstance(components[comp_name], Mass_connector):
    #             ax.plot(x, y, 'o', markersize=12, color='black', fillstyle='none', markeredgecolor='white')
    #         else:
    #             ax.text(x, y, comp_name, ha='center', va='center', fontsize=8, bbox=dict(facecolor='white', edgecolor='black'))
    
    #     plt.show()
        
#%%

"""
See for add source_sink and the connection
"""

#%% Create subcycles

wf = "Cyclopentane"
h_fluid = "INCOMP::T66"
c_fluid = "Water"

ORC = Thermal_System(wf,h_fluid,c_fluid)

#%% Create and add components

pump = Pump_cst_eff_is()
ORC.add_component(pump, "Pump",1,1,[wf])

evap_geom = Plate_HX_Geom_SWEP()
evap_geom.set_parameters_SWEP("B20Hx24/1P")

evap = Plate_HeatExchanger()
ORC.add_component(evap, "Evaporator",2,2,[wf,h_fluid])

EXP = Expander_SE()
ORC.add_component(EXP, "Expander",1,1,[wf])

cond_geom = Plate_HX_Geom_SWEP()
cond_geom.set_parameters_SWEP("B20Hx24/1P")

cond = Plate_HeatExchanger()
ORC.add_component(cond, "Condenser",2,2,[wf,c_fluid])

#%% Set Component Parameters

# Pump
pump.set_parameters(epsilon_is = 0.6, epsilon_vol = 0.8, V_s = 1e-6, V = 1.4e-3, N_pp = 50)

# Expander
EXP.set_parameters(**{
        'AU_amb': 0.6740,
        'AU_su_n': 50.0336,
        'AU_ex_n': 94.0170,
        'd_su1': 0.0097,
        'm_dot_n': 0.075,
        'A_leak': 2.6894e-06,
        'W_dot_loss_0': 10,
        'alpha': 1.2537e-05,
        'rv_in': 2.2,
        'V_s': 3.895454545454545e-05,
        'C_loss': 7.9529e-07,
        'T_amb' : 273.15 +15,
        'N_rot' : 5000
    })

# Evaporator
evap.set_parameters(**{
        'geom': evap_geom,
        'n_disc': 30,
        'wf_T': "C",
        'flow_type': "Counter_flow",
        'htc_type': "Correlation",
        'H_DP_ON': True,
        'C_DP_ON': True,
    })

# HX.set_parameters(point_su = [C_su,H_su], geom = HX_geom, n_disc = n_disc, wf_T = wf_T, flow_type = flow_type, htc_type = htc_type, H_DP_ON = DP_H_ON, C_DP_ON = DP_C_ON)

# Condenser
cond.set_parameters(**{
        'geom': cond_geom,
        'n_disc': 30,
        'wf_T': "H",
        'flow_type': "Counter_flow",
        'htc_type': "Correlation",
        'H_DP_ON': True,
        'C_DP_ON': True,
    })

#%% Link components

ORC.link_components("Pump", 1, "Evaporator", 1, wf) # Connect output 1 of pump to the input 1 of evaporator
ORC.link_components("Evaporator", 1, "Expander", 1,wf)
ORC.link_components("Expander", 1, "Condenser", 1,wf)
ORC.link_components("Condenser", 1, "Pump", 1,wf)

#%% Create and add sources/sinks

oil_source = 0
oil_sink = 0

water_source = 0
water_sink = 0

ORC.add_source(oil_source,"Oil_source_1",h_fluid,'Evaporator',2)
ORC.add_sink(oil_sink,"Oil_sink_1",h_fluid,'Evaporator',2)
ORC.add_source(water_source,"Water_source_1",c_fluid,'Condenser',2)
ORC.add_sink(water_sink,"Water_sink_1",c_fluid,'Condenser',2)

#%% Check the cycle connections
ORC.check_connections()

#%% Set known parameters : Cyclopentane loop

# After Evaporator
P_high = 36e5
ORC.set_ex(wf, 'Evaporator', 0, 'T', 273.15 + 240)
ORC.set_ex(wf, 'Evaporator', 0, 'P', P_high)
ORC.set_ex(wf, 'Evaporator', 0, 'M_dot', 0.014)

# After Turbine
P_low = 1e5
ORC.set_ex(wf, 'Expander', 0, 'P', P_low)

# After Condenser
DT_subcool = 5 # [K]
T_sat_cond = PropsSI('T','P',P_low,'Q',1,'Cyclopentane')

ORC.set_ex(wf, 'Condenser', 0, 'T', T_sat_cond-DT_subcool)

# After Pump
ORC.set_ex(wf, 'Pump', 0, 'P', P_high)

#%% Set known parameters : Water loop

# Before Condenser
ORC.set_su(c_fluid, 'Condenser', 1, 'T', 273.15 + 10)
ORC.set_su(c_fluid, 'Condenser', 1, 'P', 2e5)
ORC.set_su(c_fluid, 'Condenser', 1, 'M_dot', 0.4)

#%% Set known parameters : Oil loop

# Before Evaporator
ORC.set_su(h_fluid, 'Evaporator', 1, 'T', 273.15 + 245)
ORC.set_su(h_fluid, 'Evaporator', 1, 'P', 2e5)
ORC.set_su(h_fluid, 'Evaporator', 1, 'M_dot', 0.2)

# ORC.set_su('Cyclopentane', 'Pump', 0, 'P', 0.9*1e5)
# ORC.set_su('Cyclopentane', 'Pump', 0, 'M_dot', 0.014)

#%% See if the components are well parametrized
ORC.check_parametrized()

#%% Solve the cycle
ORC.solve()

#%% Show Results
ORC.show_cycle_info(wf)

#%% Plot the cycle

# ORC.plot_cycle(wf)  
