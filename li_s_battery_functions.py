# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 08:26:06 2019

@author: dkorff
"""

import cantera as ct
import numpy as np
from math import pi
from li_s_battery_inputs import inputs
from li_s_battery_init import cathode as cat

def dst(s1, s2, D_eff, dyInv):
    F = ct.faraday; R = ct.gas_constant; T = inputs.T
    
    C_0 = (s1['C_tot'] + s2['C_tot'])*0.5
    C_k = (s1['C_k'] + s2['C_k'])*0.5
    z_k = inputs.z_k_el
    
    N_io = np.zeros_like(s1['C_k'])
    N_io = (-D_eff*(s2['C_k'] - s1['C_k'])*dyInv
            -D_eff*C_k*(z_k*F/R/T)*(s2['phi_el'] - s1['phi_el'])*dyInv)
                
    i_io = np.dot(N_io, z_k)*F
    
    return N_io, i_io

"""========================================================================="""

def set_state(SV, offset, ptr):
    
    state = {}
    
    state['phi_ed'] = SV[offset + ptr['phi_ed']]
    state['phi_dl'] = SV[offset + ptr['phi_dl']]
    state['phi_el'] = SV[offset + ptr['phi_ed']] - SV[offset + ptr['phi_dl']]
    state['C_tot'] = sum(SV[offset + ptr['rho_k_el']])
    state['C_k'] = SV[offset + ptr['rho_k_el']]
    state['X_k'] = SV[offset + ptr['rho_k_el']]/sum(SV[offset + ptr['rho_k_el']])
    
#    state['np_S'] = SV[offset + ptr['np_S8']]
#    state['np_L'] = SV[offset + ptr['np_Li2S']]
#    state['eps_S8'] = max(SV[offset + ptr['eps_S8']], 1e-25)
#    state['eps_Li2S'] = max(SV[offset + ptr['eps_Li2S']], 1e-25)
#    
#    state['A_S'] = 3*state['eps_S8']**(2/3)/(3*cat.V_0/2/pi/state['np_S'])**(1/3)
#    state['A_L'] = 3*state['eps_Li2S']/(3*state['eps_Li2S']*cat.V_0/2/pi/state['np_L'])**(1/3)
#    
#    state['r_S'] = 3*state['eps_S8']/state['A_S']
#    state['r_L'] = 3*state['eps_Li2S']/state['A_L']
#    
#    state['A_C'] = inputs.A_C_0 \
#                 - (pi*state['np_S']*state['r_S']**2)/cat.V_0 \
#                 - (pi*state['np_L']*state['r_L']**2)/cat.V_0
    
    return state

"""========================================================================="""

def set_state_sep(SV, offset, ptr):
    
    state = {}
    
    state['phi_el'] = SV[offset + ptr['phi']]
    state['C_tot'] = sum(SV[offset + ptr['rho_k_el']])
    state['C_k'] = SV[offset + ptr['rho_k_el']]
    state['X_k'] = SV[offset + ptr['rho_k_el']]/sum(SV[offset + ptr['rho_k_el']])
    
    return state