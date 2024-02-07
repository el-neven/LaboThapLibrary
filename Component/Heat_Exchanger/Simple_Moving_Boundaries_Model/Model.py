import CoolProp.CoolProp as cp
from CoolProp.CoolProp import PropsSI
import math
import numpy as np
from scipy.optimize import least_squares


def U_DittusBoelter(m, heat_transfer_direction, D, fluid, P, T):
    if heat_transfer_direction == 'heated':
        a = 0.4
    elif heat_transfer_direction == 'cooled':
        a = 0.3
    
    if len(locals()) > 6:
        # case 1 phase heated or cooled
        mu = PropsSI('V', 'T', T, 'P', P, fluid)
        cp = PropsSI('CPMASS', 'T', T, 'P', P, fluid)
        # k = PropsSI('L', 'T', T, 'P', P, fluid)
        k=0.08674 

        Re = m * 4 / (math.pi * mu * D)
        Pr = mu * cp / k
        U = 0.023 * k * Re**0.8 * Pr**a / D
    elif len(locals()) == 6:
        # case 2 phases
        mu = (cp.PropsSI('V', 'P', P, 'Q', 0, fluid) + cp.PropsSI('V', 'P', P, 'Q', 1, fluid)) / 2
        cp = (cp.PropsSI('C', 'P', P, 'Q', 0, fluid) + cp.PropsSI('C', 'P', P, 'Q', 1, fluid)) / 2
        k = (cp.PropsSI('L', 'P', P, 'Q', 0, fluid) + cp.PropsSI('L', 'P', P, 'Q', 1, fluid)) / 2
        Re = m * 4 / (math.pi * mu * D)
        Pr = mu * cp / k
        U = 0.023 * k * Re**0.8 * Pr**a / D
    
    return U


def U_Thonon(m, D, fluid, P, T):
    mu = PropsSI('V', 'T', T, 'P', P, fluid)
    cp = PropsSI('CPMASS', 'T', T, 'P', P, fluid)
    k=0.08674 
    # k = PropsSI('L', 'T', T, 'P', P, fluid)
    Re = m * 4 / (math.pi * mu * D)
    Pr = mu * cp / k
    U = 0.2946 * k * Re**0.7 * Pr**(1/3) / D
    return U


def U_Gnielinski_calibrated(m, D, fluid, P):
    T = PropsSI('T', 'P', P, 'Q', 0, fluid)
    try:
        mu = PropsSI('V', 'T', T, 'P', P, fluid)
    except:
        mu = (PropsSI('V', 'P', P, 'Q', 0, fluid) + PropsSI('V', 'P', P, 'Q', 1, fluid)) / 2
    
    try:
        cp = PropsSI('CPMASS', 'T', T, 'P', P, fluid)
    except:
        cp = (PropsSI('C', 'P', P, 'Q', 0, fluid) + PropsSI('C', 'P', P, 'Q', 1, fluid)) / 2

    # k = PropsSI('L', 'T', T, 'P', P, fluid)
    k=0.08674 
    Re = m * 4 / (math.pi * mu * D)
    Pr = mu * cp / k
    f = (0.79 * math.log(Re) - 1.64) ** (-2)
    c = 4.616163048309070
    Nu = c * (f / 8 * (Re - 1000) * Pr) / (1 + 12.7 * (f / 8) ** 0.5 * (Pr ** (2 / 3) - 1))
    U = k * Nu / D
    return U

def U_Cooper_calibrated(Q, HX_A, P, fluid):

    pr = P / cp.PropsSI('Pcrit', fluid)
    M = cp.PropsSI('M', 'C', 0, ' ', 0, fluid)
    c = [0.978552748683140, 1.07700234277466]
    U = c[0] * (Q / HX_A) ** (0.67 * c[1]) * 55 * pr ** (0.12 - 0.2 * math.log(pr)) * (-math.log(pr)) ** (-0.55 * c[1]) * M ** (-0.5)
    return U


class HX_model_HP_CB():

    def __init__(self): #   dT_sub_or_sup, fluid_wf, P_sf, T_sf_in, glide_sf, fluid_sf, x_eva_in):
        
        self.calculable = False
        self.parametrized = False
        self.defined = False

        "Initialize connectors"
        self.point_su1 = None
        self.point_su2 = None
        self.point_ex1 = None
        self.point_ex2 = None

        "Parameters"
        self.HX_type = None
        self.HX_D = None
        self.HX_A = None
        self.glide_sf = None
        #self.dT_sub_or_sup = None
        self.min_pinch = None
        self.max_pinch = None

    def update_inputs(self, point_su1, point_su2, point_ex1, point_ex2):

        self.point_su1 = point_su1
        self.point_su2 = point_su2
        self.point_ex1 = point_ex1
        self.point_ex2 = point_ex2

        "Inputs"
        self.fluid_wf = point_su1.fluid
        self.T_su_wf = point_su1.T
        # self.P_su_wf = point_su1.p
        self.m_dot_wf = point_su1.m_dot

        self.fluid_sf = point_su2.fluid
        self.T_su_sf = point_su2.T
        self.P_su_sf = point_su2.p
        self.m_dot_sf = point_su2.m_dot

        self.T_ex_wf = point_ex1.T

        # self.point_su1 = point_su1
        # self.point_su2 = point_su2

        # self.point_ex1 = point_ex1
        # self.point_ex2 = point_ex2

        # self.Required_inputs = [self.point_su1.T, self.T_ex_wf, self.m_dot_wf, self.T_su_sf, self.P_su_sf, self.m_dot_sf]

        self.check_calculable()

    def check_calculable(self):
        self.Required_inputs = [self.point_su1.T, self.point_ex1.T, self.point_su1.p, self.point_su1.m_dot, self.point_su2.T, self.point_su2.p, self.point_su2.m_dot, self.point_su1.fluid]

        if all(Input is not None for Input in self.Required_inputs):
            self.calculable = True

    def check_parametrized(self):
        if self.glide_sf != None and self.HX_type != None and self.HX_D != None and self.HX_A != None and self.glide_sf != None:
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
        # Resolution of the problem - iterative variable P_wf: saturation pressure of the working fluid

        if self.HX_type == 'evaporator':
            lb = PropsSI('P', 'T', self.T_su_sf - 2 - self.max_pinch, 'Q', 0, self.fluid_wf) #A changer ici!!
            ub = PropsSI('P', 'T', self.T_su_sf - 2 - self.min_pinch, 'Q', 0, self.fluid_wf) #A changer ici!!
            # lb = self.T_su_sf - self.T_su_wf - self.glide_sf - self.max_pinch
            # ub = self.T_su_sf - self.T_su_wf - self.glide_sf - self.min_pinch
            #sf_plot_color = 'b'
        elif self.HX_type == 'condenser':
            #print(self.T_su_sf)
            # print(self.min_pinch)
            # print(self.fluid_wf)
            lb = PropsSI('P', 'T', self.T_su_sf + 2 + self.min_pinch, 'Q', 0, self.fluid_wf)
            ub = PropsSI('P', 'T', self.T_su_sf+ 2 + self.max_pinch, 'Q', 0, self.fluid_wf)
            # lb = self.T_su_wf - self.T_su_sf - self.max_pinch
            # ub = self.T_su_wf - self.T_su_sf - self.min_pinch
            # sf_plot_color = 'r'
        
        # # Minimization problem
        # if abs(ub - lb) < 10 or (ub < lb):
        #     P_wf0 = (lb + ub) / 2
        # else:
        #     P_wf0 = np.linspace(lb, ub, 10)
        P_wf0 = (lb + ub) / 2
        #print(P_wf0)
        least_squares(self.equations, P_wf0, bounds=(lb, ub))

        print(self.P_sat_wf)
        self.point_ex1.set_fluid(self.point_su1.fluid)
        self.point_ex1.set_m_dot(self.point_su1.m_dot)
        self.point_ex1.set_T(self.T_ex_wf)
        self.point_ex1.set_p(self.P_sat_wf) #prendre en compte le pressure drops ici
        self.defined = True
        # for ii in range(len(P_wf0)):
        #     self.results_temp = self.equations(P_wf0[ii])

        # #Pas encore vérifier jusqu'ici 
        #     self.flag_tau[ii] = results_temp['flag_tau']
        #     self.res_temp[ii] = results_temp['res']
        
        # self.res_temp[self.flag_tau == 1] = []
        # if len(self.res_temp) == 0:   # Check compatibility of boundary conditions
        #     self.results['P_wf'] = float('nan')
        #     self.results['flag_check'] = 0
        # else:
        #     P_wf0[self.flag_tau == 1] = []
        #     min_indx = self.res_temp.index(min(self.res_temp))
        #     self.results = self.moving_boundaries_HX_model(P_wf0[min_indx])
        #     self.results['P_wf'] = P_wf0[min_indx]
        #     self.results['flag_check'] = 1

        
    

    def equations(self, P_sat_wf):
        
        "Inputs"
        fluid_wf = self.point_su1.fluid
        T_su_wf = self.point_su1.T
        #print(T_su_wf-273.15)
        # P_sat_wf = self.point_su1.p
        # x_su_wf = self.point_su1.x
        m_dot_wf = self.point_su1.m_dot

        T_ex_wf = self.point_ex1.T

        fluid_sf = self.point_su2.fluid
        T_su_sf = self.point_su2.T
        P_su_sf = self.point_su2.p
        m_dot_sf = self.point_su2.m_dot
        
        # TEMPERATURES AND THERMAL POWERS
        # evaporating/consensing phase - wf side
        T_wf_sat = PropsSI('T', 'P', P_sat_wf, 'Q', 0.5, fluid_wf)
        #print(T_su_wf, T_wf_sat)
        h_wf_sat_vap = PropsSI('H', 'P', P_sat_wf, 'Q', 1, self.fluid_wf)
        x_su_wf = PropsSI('Q', 'P', P_sat_wf, 'T', T_su_wf, self.fluid_wf)

        #If two-phase, inlet temperature is the saturation temperature
        if x_su_wf >= 0 and x_su_wf < 1:
            #Enthalpy of the saturated liquid
            h_wf_sat_liq = PropsSI('H', 'P', P_sat_wf, 'Q', self.x_su_wf, self.fluid_wf)
            T_su_wf = T_wf_sat
        else:
            h_wf_sat_liq = PropsSI('H', 'P', P_sat_wf, 'Q', 0, self.fluid_wf)
    
        # case depending temperatures
        if self.HX_type == 'evaporator':
            sf_ht_dir = 'cooled'
            wf_ht_dir = 'heated'
            
            #In the case of an evaporatore, the secondary fluid is cooled
            T_sf_h = T_su_sf
            T_sf_c = T_sf_h - self.glide_sf

            T_wf_h = T_ex_wf # Est-ce qu'on peut se débarrasser de ça???
            T_wf_c = T_su_wf
            # self.T_ex_wf = T_wf_h

        elif self.HX_type == 'condenser':
            sf_ht_dir = 'heated'
            wf_ht_dir = 'cooled'
            
            #In the case of a condenser, the secondary fluid is heated
            T_sf_c = T_su_sf
            T_sf_h = T_sf_c + self.glide_sf

            T_wf_c = T_ex_wf # Est-ce qu'on peut se débarrasser de ça???
            T_wf_h = T_su_wf
            # self.T_ex_wf = T_wf_c
    
        # Energy balance - wf side
        h_wf_h = PropsSI('H', 'T', T_wf_h, 'P', P_sat_wf, fluid_wf) #Enthalpy on the wf hot side
        h_wf_c = PropsSI('H', 'T', T_wf_c, 'P', P_sat_wf, fluid_wf) #Enthalpy on the wf cold side
        Q_sub = max(0, m_dot_wf * (h_wf_sat_liq - h_wf_c)) #Heat transfered in the subcooled zone
        Q_sat = max(0, m_dot_wf * (h_wf_sat_vap - h_wf_sat_liq)) #Heat transfered in the two-phase zone
        Q_sup = max(0, m_dot_wf * (h_wf_h - h_wf_sat_vap)) #Heat transfered in the superheated zone
        self.Q = Q_sub + Q_sat + Q_sup #Sum of the three zones
        # print(h_wf_sat_liq, h_wf_c)
        #print(T_wf_sat, T_wf_c)
        # Energy balance - sf side
        T_sf_mean = (T_sf_h + T_sf_c) / 2 #Mean temperature on the sf side
        cp_sf = PropsSI('C', 'T', T_sf_mean, 'P', P_su_sf, fluid_sf)
        self.m_dot_sf = self.Q / (cp_sf * (T_sf_h - T_sf_c)) #The code gives you the mass flow rate of the sf
        T_sf_ch = T_sf_c + Q_sub / (self.m_dot_sf * cp_sf) #Temperature of the sf at the end of the subcooled zone
        T_sf_hc = T_sf_h - Q_sup / (self.m_dot_sf * cp_sf) #Temperature of the sf at the end of the superheated zone
        # print(T_sf_ch, T_sf_c, Q_sub)
        # FLUID MASS DISTRIBUTION
        # Heat tranfer coefficients 1P - sf side: Dittus_Boelter
        U_sf = U_DittusBoelter(self.m_dot_sf, sf_ht_dir, self.HX_D, fluid_sf, P_su_sf, T_sf_mean)
        
        # Heat tranfer coefficients - wf
        if self.HX_type == 'condenser':
            self.U_wf_h = U_Thonon(m_dot_wf, self.HX_D, fluid_wf, P_sat_wf, (T_wf_h + T_wf_sat) / 2) #Heat transfer coefficient at the beginning of the superheated zone
            self.U_wf_c = U_Thonon(m_dot_wf, self.HX_D, fluid_wf, P_sat_wf, (T_wf_c + T_wf_sat) / 2) #Heat transfer coefficient at the end of the subcooled zone
            self.U_wf_2P = U_Cooper_calibrated(self.Q, self.HX_A, P_sat_wf, fluid_wf) #Heat transfer coefficient in the two-phase zone
        elif self.HX_type == 'evaporator':
            self.U_wf_h = U_Thonon(m_dot_wf, self.HX_D, fluid_wf, P_sat_wf, (T_wf_h + T_wf_sat) / 2) #Heat transfer coefficient at the beginning of the superheated zone
            self.U_wf_c = U_Thonon(m_dot_wf, self.HX_D, fluid_wf, P_sat_wf, (T_wf_c + T_wf_sat) / 2) #Heat transfer coefficient at the end of the subcooled zone
            self.U_wf_2P = U_Gnielinski_calibrated(m_dot_wf, self.HX_D, fluid_wf, P_sat_wf) #Heat transfer coefficient in the two-phase zone
        
        self.U_sub = (1 / U_sf + 1 / self.U_wf_c) ** (-1) #Heat transfer coefficient in the subcooled zone (addition of the two heat transfer coefficients in series)
        self.U_sat = (1 / U_sf + 1 / self.U_wf_2P) ** (-1) #Heat transfer coefficient in the two-phase zone (addition of the two heat transfer coefficients in series)
        self.U_sup = (1 / U_sf + 1 / self.U_wf_h) ** (-1) #Heat transfer coefficient in the superheated zone (addition of the two heat transfer coefficients in series)
    
        # Mean logarithm temperature difference
        if self.HX_type == 'evaporator':
            tau_sub = [(T_sf_c - T_wf_c), (T_sf_ch - T_wf_sat)] #Temperature difference between the secondary fluid and the working fluid in the subcooled zone
            tau_sat = [(T_sf_ch - T_wf_sat), (T_sf_hc - T_wf_sat)] #Temperature difference in the two-phase zone
            tau_sup = [(T_sf_hc - T_wf_sat), (T_sf_hc - T_wf_h)] #Temperature difference in the superheated zone
        elif self.HX_type == 'condenser':
            tau_sub = [-(T_sf_c - T_wf_c), -(T_sf_ch - T_wf_sat)] #Temperature difference between the secondary fluid and the working fluid in the subcooled zone
            tau_sat = [-(T_sf_ch - T_wf_sat), -(T_sf_hc - T_wf_sat)] #Temperature difference in the two-phase zone
            tau_sup = [-(T_sf_hc - T_wf_sat), -(T_sf_hc - T_wf_h)] #Temperature difference in the superheated zone

        if (T_sf_c - T_wf_c)<0 or (T_sf_ch - T_wf_sat)<0 or (T_sf_hc - T_wf_sat)<0:
            self.flag_tau = True
            dT_ml_sub = 1
            dT_ml_sat = 1
            dT_ml_sup = 1
            #print('coucou')
        else:
            #print('quoi?')
            dT_ml_sub = (max(tau_sub) - min(tau_sub)) / np.log(max(tau_sub) / min(tau_sub)) #Mean logarithm temperature difference in the subcooled zone
            dT_ml_sat = (max(tau_sat) - min(tau_sat)) / np.log(max(tau_sat) / min(tau_sat)) #Mean logarithm temperature difference in the two-phase zone
            dT_ml_sup = (max(tau_sup) - min(tau_sup)) / np.log(max(tau_sup) / min(tau_sup)) #Mean logarithm temperature difference in the superheated zone
            # print(dT_ml_sub, dT_ml_sat, dT_ml_sup)
        self.tau = tau_sub + tau_sat + tau_sup #Sum of the three zones
        # self.flag_tau = any(t < -0.5 for t in self.tau)
    
        # Energy balance UA
        self.A_sub = Q_sub / (self.U_sub * dT_ml_sub) #Heat transfer surface in the subcooled zone
        self.A_sup = Q_sup / (self.U_sup * dT_ml_sup) #Heat transfer surface in the superheated zone
        self.A_sat = abs(self.HX_A - self.A_sub - self.A_sup) #Heat transfer surface in the two-phase zone
        
        # Closure equations
        Q_sat_calc = self.A_sat * self.U_sat * dT_ml_sat #Heat transfered in the two-phase zone
        self.res = abs(Q_sat - Q_sat_calc)
        self.P_sat_wf = P_sat_wf

        #print(self.res)
        return self.res
        # Results
        self.point_ex1.set_fluid(self.point_su1.fluid)
        self.point_ex1.set_m_dot(self.point_su1.m_dot)
        self.point_ex1.set_T(T_ex_wf)
        self.point_ex1.set_p(self.point_su1.p)
        print(self.point_ex1.h)
        self.defined = True
        
        # # Results for TQ plot
        # out['plot'] = {}
        # out['plot']['Q_wf'] = [0, Q_sub, Q_sub + Q_sat, out['Q']]
        # out['plot']['Q_sf'] = [0, out['Q']]
        # out['plot']['T_wf'] = [T_wf_c, T_wf_sat, T_wf_sat, T_wf_h]
        # out['plot']['T_sf'] = [T_sf_c, T_sf_h]

    # def heat_transfer_diagram(self, results):
    #     if self.TQ_flag == 1:
    #         HEAT_TRANSFER_DIAGRAM = figure()
    #         Size = 18
    #         font = 'Arial'
    #         HEAT_TRANSFER_DIAGRAM.set_position([100, 100, 600, 420])
    #         plot_ax = HEAT_TRANSFER_DIAGRAM.add_subplot(111)
    #         plot_ax.hold(True)
    #         plot_ax.set_ylabel('Temperature (°C)', fontname=font, fontsize=Size, fontweight='bold')
    #         plot_ax.set_xlabel('Termal power (kW)', fontname=font, fontsize=Size, fontweight='bold')
    #         plot_ax.set_title(param['HX_type'] + ' heat transfer diagram')
    #         plot_ax.box(True)
    #         plot_ax.set_fontname('Arial')
    #         plot_ax.set_fontsize(Size - 2)
    #         plot_ax.grid(True, linestyle=':', color=[0, 0, 0], alpha=0.5)
    #         plot_ax.plot(results['plot']['Q_wf'] / 1e3, results['plot']['T_wf'] - 273.15, 'g', linewidth=2)
    #         plot_ax.plot(results['plot']['Q_sf'] / 1e3, results['plot']['T_sf'] - 273.15, sf_plot_color, linewidth=2)
    #         plot_ax.text(np.mean(results['plot']['Q_wf'] / 1e3), results['plot']['T_wf'][2] - 273.15, str(round(results['P_wf'] / 100)) + ' bar', backgroundcolor='w')
    #         plot_ax.text(np.mean(results['plot']['Q_sf'] / 1e3), np.mean(results['plot']['T_sf']) - 273.15, str(round(results['m_dot_sf'], 1)) + ' l/s', backgroundcolor='w')
    #         plot_ax.hold(False)
        
    #     return results
    

