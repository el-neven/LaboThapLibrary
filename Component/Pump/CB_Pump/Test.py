from Model import PumpModel
import sys
sys.path.append('C:/Users/Elise/OneDrive - Universite de Liege/Documenten/Thesis/LaboThapLibrary')

from Port.Mass_connector import Mass_connector

"Define ports"
point_su = Mass_connector()
point_ex = Mass_connector()

point_su.set_fluid('R1234yf')
point_su.set_m_dot(0.01)
point_su.set_T(230.92)
point_su.set_p(200000)

point_ex.set_fluid('R1234yf')
point_ex.set_p(400000)


m_dot = 1.2
Head = 50
Pump = PumpModel()
Pump.set_inputs(point_su, point_ex)
Pump.solve()
