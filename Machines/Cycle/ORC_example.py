# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 15:23:02 2024

@author: Elise
"""
import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/Python_Library')

# Modèle des composants principaux, sans perte de charge et perte à l'ambiance dans les tuyères et composants auxilliaires

from CoolProp.CoolProp import PropsSI

# import sys
# sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/Python_Library/Component/Heat_Exchanger/Plate_Heat_Exchanger/Moving_Boundary_Model/Modules')

from Port.Mass_connector import Mass_connector
from Port.Heat_connector import Heat_connector
from Component.Expander.Semi_Empirical_Model_VL.Simulation_Model import Expander_SE
from Component.Heat_Exchanger.Simple_Moving_Boundaries_Model.Model import HX_model_HP_CB
from Component.Pump.Constant_Efficiency.Model import Pump_cst_eff_is
from Cycle import Cycle
from scipy.optimize import fsolve
from scipy.optimize import root
from scipy.optimize import minimize
from scipy.optimize import least_squares

"""
This script is used to test the Cycle class. The cycle is defined by the components and the number of points. The test is carried out with the ORC example of the POC1.

The inputs are the evaporation and condensation pressure guesses and the bounds of the optimization problem.
No pressure losses are taken into account in the cycle.
The models used are: semi-empirical model for the expander, simple moving boundaries model for the heat exchangers and constant efficiency model for the pump.
"""

def Cycle_resolution(P_ev_guess, P_cd_guess, bounds): #Function with the optimization problem


    "1. Define the componants of the cycle and their parameters"
    list_of_components = []

    # Point 1->2: Expander
    EXP = Expander_SE()

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
        'C_loss': 7.9529e-07
    })

    list_of_components = list_of_components + [EXP]

    #Point 2->3: Condenser
    COND = HX_model_HP_CB()

    COND.set_parameters(**{
        'HX_type': 'condenser',
        'HX_D': 0.02,
        'HX_A': 2.38,
        'min_pinch': 0.5,
        'max_pinch': 20,
        'glide_sf': 10,
    })

    list_of_components = list_of_components + [COND]

    #Point 3->4: Pump
    PUMP = Pump_cst_eff_is()
    PUMP.set_parameters(epsilon_is = 0.5, epsilon_vol = 0.8, V_s = 1e-6, V = 1.4e-3)

    list_of_components = list_of_components + [PUMP]

    #Point 4->1: Evaporator
    EVAP = HX_model_HP_CB()

    EVAP.set_parameters(**{
        'HX_type': 'evaporator',
        'HX_D': 0.02,
        'HX_A': 2.38,
        'min_pinch': 0.5,
        'max_pinch': 20,
        'glide_sf': 10,
    })

    list_of_components = list_of_components + [EVAP]

    "2. Cycle definition"

    cycle = Cycle(4, list_of_components) #The cycle is defined by the number of points and the list of components

    "3. Inputs of the model"

    # "Inputs"
    Delta_T_sh = 5 # [K]
    Delta_T_sc = 5 # [K]

    "Work connectors"
    N_exp = 5000
    N_pp = 200

    "Heat connectors"
    Amb = Heat_connector()

    T_amb = 293
    Amb.set_T1(T_amb)

    # Define the system of equations
    def equations(P):

        P_ev, P_cd = P
        "4. Set the variable defined to False so that the cycle can be solved"
        for component in list_of_components:
            component.defined = False
        
        "5. Set the values of the connectors"
        # Point 0
        cycle.point[0].set_p(P_ev) #[Pa]
        T_0 = PropsSI('T','P',P_ev,'Q',0.5,'R1233zd(E)')+Delta_T_sh
        cycle.point[0].set_T(T_0) #[K]
        cycle.point[0].set_fluid('R1233zd(E)')
        h_0_init = cycle.point[0].h

        # Point 1
        cycle.point[1].set_p(P_cd) #[Pa]
        cycle.point[1].set_fluid('R1233zd(E)') #[Pa]

        #Point 2
        cycle.point[2].set_p(P_cd) #[Pa]
        T_2 = PropsSI('T','P',P_cd,'Q',0.5,'R1233zd(E)')-Delta_T_sc
        cycle.point[2].set_T(T_2) #[K]
        cycle.point[2].set_fluid('R1233zd(E)')
        h_2_init = cycle.point[2].h

        #Point 4
        cycle.point[3].set_p(P_ev) #[Pa]
        cycle.point[3].set_fluid('R1233zd(E)') #[Pa]

        "Heat source cycle"
        htf_su = Mass_connector()
        htf_su.set_p(2e5) #[Pa]
        htf_su.set_T(72.5367743441399+273.15) #[K]
        htf_su.set_m_dot(0.199615970833042) #[kg/s]
        htf_su.set_fluid('Water')

        htf_ex = Mass_connector()

        "Cold source cycle"
        ctf_su = Mass_connector()
        ctf_su.set_p(1e5) #[Pa]
        ctf_su.set_T(14.3799641808547+273.15) #[J/kg]
        ctf_su.set_m_dot(0.267194284887110) #[kg/s]
        ctf_su.set_fluid('Water')

        ctf_ex = Mass_connector()

        #Point 1->2: Expander
        EXP.update_inputs(cycle.point[0], cycle.point[1], N_exp, Amb.T_1)

        #Point 2->3: Condenser
        COND.update_inputs(cycle.point[1], ctf_su, cycle.point[2], ctf_ex)

        #Point 3->4: Pump
        PUMP.update_inputs(cycle.point[2], cycle.point[3], N_pp)

        #Point 4->1: Evaporator
        EVAP.update_inputs(cycle.point[3], htf_su, cycle.point[0], htf_ex)

        "5. Solver le cycle"
        cycle.solve()

        "Residues"
        res_h_0 = abs((cycle.point[0].h - h_0_init)/h_0_init)
        res_h_2 = abs((cycle.point[2].h - h_2_init)/h_2_init)
        res = [res_h_0, res_h_2]
        print(res)
        return res
    
    P_guess = [P_ev_guess, P_cd_guess]

    result = least_squares(equations, P_guess, bounds=bounds)

    P_solution = result.x
    P_ev_solution, P_cd_solution = P_solution

    # Return the optimized values
    return P_ev_solution, P_cd_solution


# Set the guess and boundary constraints
T_su_sf = 73+273.15

P_ev_lb = PropsSI('P', 'T', T_su_sf - 5 - 20, 'Q', 0, 'R245fa') #A changer ici!!
P_ev_ub = PropsSI('P', 'T', T_su_sf - 5 - 0.5, 'Q', 0, 'R245fa') #A changer ici!!

P_cd_lb = PropsSI('P', 'T', T_su_sf + 5 + 0.5, 'Q', 0, 'R245fa')
P_cd_ub = PropsSI('P', 'T', T_su_sf + 5 + 20, 'Q', 0, 'R245fa')

P_ev_guess = (P_ev_lb+P_ev_ub)/2
P_cd_guess = (P_cd_lb+P_cd_ub)/2

bounds = [(P_ev_lb, P_cd_lb), (P_ev_ub, P_cd_ub)]

P_ev_optimized, P_cd_optimized = Cycle_resolution(P_ev_guess, P_cd_guess, bounds)

