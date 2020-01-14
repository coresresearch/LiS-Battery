# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 08:26:06 2019

@author: dkorff
"""

import cantera as ct
import numpy as np
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
    
#    N_io = N_io*np.array((1., 1., 1., 1., 0., 0., 0., 0., 0., 0.))
        
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
    
    return state

"""========================================================================="""

def set_state_sep(SV, offset, ptr):
    
    state = {}
    
    state['phi_el'] = SV[offset + ptr['phi']]
    state['C_tot'] = sum(SV[offset + ptr['rho_k_el']])
    state['C_k'] = SV[offset + ptr['rho_k_el']]
    state['X_k'] = SV[offset + ptr['rho_k_el']]/sum(SV[offset + ptr['rho_k_el']])
    
    return state