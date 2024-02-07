# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 14:39:37 2023

@author: Elise
"""

import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/Python_Library')


from Port.Mass_connector import Mass_connector
from Simulation_Model import Expander_SE

import numpy as np
import matplotlib.pyplot as plt

"Define ports"
point_su = Mass_connector()
point_ex = Mass_connector()

point_su.set_fluid('R1234yf')
#point_su.set_m_dot(0.01)
point_su.set_T(288.92)
point_su.set_p(400000)

point_ex.set_p(400000/2)

"Work connector (pas encore coder!!)"
N_exp = 5000/1.4

"Heat connector (pas encore coder!!)"
T_amb = 293
"Define class"
EX = Expander_SE(point_su, point_ex, N_exp, T_amb)

EX.set_parameters(**{
    'AU_amb': 5,
    'AU_su_n': 14,
    'AU_ex_n': 20,
    'd_su1': 8e-3,
    'm_dot_n': 0.09,
    'A_leak': 4e-7,
    'W_dot_loss_0': 10,
    'alpha': 0.01,
    'rv_in': 2,
    'V_s': 0.000034,
    'C_loss': 0
})

EX.solve()
# #---------------------------------------------------------------------------------
# # Graph varying the pressure ratio

# # Create lists to store results
# rp_values = np.linspace(1.3, 6, 25)
# epsilon_is = []

# # Loop through different values of CP.set_rp
# for rp in rp_values:
    
#     EX.set_rp(rp)
    
#     EX.solve()
#     epsilon_is.append(EX.epsilon_is)


# plt.plot(rp_values, epsilon_is, linewidth=1.9)

# plt.xlabel(r'$\mathrm{r_{p}}$ [-]', fontsize=16)
# plt.ylabel(r'$\mathrm{\epsilon_{is}}$ [-]', fontsize=16)

# plt.ylim([0, 1.05])
# plt.xlim([1.4, 6])

# plt.grid(True)
# plt.show()




