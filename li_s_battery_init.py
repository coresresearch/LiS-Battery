# -*- coding: utf-8 -*-
"""
Created on Wed May  8 10:39:12 2019

@author: dkorff

This is the object initialization file for the Li-s model. It imports input
values from li_s_battery_inputs.py and initializes the necessary objects to
run the simulations
"""

import numpy as np
import cantera as ct
import importlib
from math import pi

#import li_s_battery_inputs
#importlib.reload(li_s_battery_inputs)
from li_s_battery_inputs import inputs


"Import cantera objects - this step is the same regardless of test type"
#anode_obj = ct.Solution(inputs.ctifile, inputs.anode_phase)
elyte_obj = ct.Solution(inputs.ctifile, inputs.elyte_phase)
sulfur_obj = ct.Solution(inputs.ctifile, inputs.cat_phase1)
Li2S_obj = ct.Solution(inputs.ctifile, inputs.cat_phase2)
carbon_obj = ct.Solution(inputs.ctifile, inputs.cat_phase3)
conductor_obj = ct.Solution(inputs.ctifile, inputs.metal_phase)
lithium_obj = ct.Solution(inputs.ctifile, inputs.an_phase)

#anode_s_obj = ct.Interface(inputs.ctifile, inputs.anode_surf_phase,
#                           [anode_obj, elyte_obj, conductor_obj])
sulfur_el_s = ct.Interface(inputs.ctifile, inputs.sulfur_elyte_phase,
                             [sulfur_obj, elyte_obj, conductor_obj])
Li2S_el_s = ct.Interface(inputs.ctifile, inputs.Li2S_elyte_phase,
                             [Li2S_obj, elyte_obj, conductor_obj])
carbon_el_s = ct.Interface(inputs.ctifile, inputs.graphite_elyte_phase,
                             [carbon_obj, elyte_obj, conductor_obj])
lithium_el_s = ct.Interface(inputs.ctifile, inputs.anode_elyte_phase,
                             [lithium_obj, elyte_obj, conductor_obj])
#Li2S_tpb = ct.Interface(inputs.ctifile, 'tpb', [Li2S_obj, Li2S_el_s, ])

if hasattr(inputs, 'C_k_el_0'):
    elyte_obj.X = inputs.C_k_el_0/np.sum(inputs.C_k_el_0)

# Set initial conditions of cantera objects
#Li2_S = inputs.Li_cat_max - inputs.SOC_0*(inputs.Li_cat_max - inputs.Li_cat_min)
#S8 = 1 - Li2_S
#X_cat_init = '{}:{}, {}:{}'.format(inputs.Li_species_cathode, Li2_S,
#                                   inputs.Vac_species_cathode, S8)
#
#cathode_obj.X = X_cat_init
#cathode_obj.TP = inputs.T, ct.one_atm
#cathode_s_obj.TP = inputs.T, ct.one_atm
#
#X_cat_0 = cathode_obj.X
#X_elyte_0 = elyte_obj.X

bc_class = getattr(inputs, 'test_type')

F = ct.faraday

#class battery():
    
#    def get_i_ext():
#        return battery.i_ext
    
#    def set_i_ext(value):
#        battery.i_ext = value
        
#    oneC = inputs.epsilon_cat*cathode_obj.density_mole*inputs.H_cat*F/3600
    
    # Calculate the actual current density. 
#    if inputs.flag_cathode == 1:
#        i_ext_amp = -inputs.C_rate*oneC_cat
       
       
"============================================================================="

class cathode():
    # Set a flag to let the solver know whether to implement this class
    flag = inputs.flag_cathode
    F = ct.faraday
    
    # Number of nodes in the y-direction
    npoints = inputs.npoints_cathode
    
    # Number of shells in the cathode particle
#        nshells = inputs.nshells_cathode
    
    # Number of state variables per node
    nVars = 2 + elyte_obj.n_species + 4
    
    # Pointers
    ptr = {}
    ptr['iFar'] = elyte_obj.species_index(inputs.Li_species_elyte)
    
    ptr['eps_S8'] = 0
    ptr['eps_Li2S'] = 1
    ptr['rho_k_el'] = 2 + np.arange(0, elyte_obj.n_species)
    ptr['phi_dl'] = ptr['rho_k_el'][-1] + 1
    ptr['phi_ed'] = ptr['rho_k_el'][-1] + 2
    ptr['np_S8'] = ptr['rho_k_el'][-1] + 3
    ptr['np_Li2S'] = ptr['rho_k_el'][-1] + 4
    
    nSV = npoints*nVars
    offsets = np.arange(0, int(nSV), int(nVars))
    
    ptr_vec = {}
    ptr_vec['eps_S8']   = ptr['eps_S8']   + offsets
    ptr_vec['eps_Li2S'] = ptr['eps_Li2S'] + offsets
    ptr_vec['rho_k_el'] = ptr['rho_k_el']
    for i in offsets[1:]:
        ptr_vec['rho_k_el'] = np.hstack((ptr_vec['rho_k_el'],i+ptr['rho_k_el']))
    ptr_vec['phi_dl']  = ptr['phi_dl'] + offsets
    ptr_vec['phi_ed']  = ptr['phi_ed']   + offsets
    ptr_vec['np_S8']   = ptr['np_S8']   + offsets
    ptr_vec['np_Li2S'] = ptr['np_Li2S'] + offsets
    
    # Store parameters as class attributes
    T = inputs.T
    C_dl = inputs.C_dl_cat
    
    # Geometric parameters
    tau = inputs.tau_cat
    r_p = inputs.r_p_cat
    d_p = inputs.d_p_cat
    dyInv = npoints/inputs.H_cat
    dy = inputs.H_cat/npoints
    H = inputs.H_cat
    V_0 = inputs.H_cat*inputs.A_cat
    
    
    if inputs.sulfur_method == 'bulk':
        m_S = inputs.m_S_0/inputs.A_cat
        m_S_0 = inputs.m_S_0
    elif inputs.sulfur_method == 'loading':
        m_S = inputs.m_S_0
        m_S_0 = inputs.m_S_0*inputs.A_cat
        
        
    omega_S = inputs.pct_w_S8_0/inputs.A_cat
    omega_C = inputs.pct_w_C_0/inputs.A_cat
    rho_S = sulfur_obj.density_mass
    rho_C = carbon_obj.density_mass
    m_solid = m_S/omega_S
    
    eps_S_0 = m_S/rho_S/H
    eps_C_0 = m_solid*omega_C/rho_C/H
    eps_L_0 = 1e-5; #A_L_0 = 1e5
    
    A_S_0 = (3*eps_S_0)/((3*eps_S_0*V_0)/(2*pi*inputs.np_S8_init))**(1/3)
    A_L_0 = (3*eps_L_0)/(3*eps_L_0*V_0/2/inputs.np_Li2S_init/pi)**(1/3)
    
#    A_S_0 = 3*(3*eps_S_0*V_0/2/pi/inputs.np_S8_init)**(-1/3)
#    A_L_0 = 3*(3*eps_L_0*V_0/2/pi/inputs.np_Li2S_init)**(-1/3)
    
    eps_el_0 = 1 - eps_S_0 - eps_C_0 - eps_L_0
    eps_pore = 1 - eps_C_0
    
#    r_C = (3*eps_C_0*inputs.A_cat*H/4/inputs.npoints_cathode/pi)**(1/3)
    r_C = 3*eps_C_0/inputs.A_C_0
#    r_C = inputs.H_cat
#    A_C_0 = 3*eps_C_0/(H/2)
    
    oneC = 16*(eps_S_0)*sulfur_obj.density_mole*H*F/3600
#    oneC = eps_S_0*H*sulfur_obj.density_mass*1675
    
    def get_i_ext():
        return cathode.i_ext
    
    def set_i_ext(value):
        cathode.i_ext = value
            
    # Calculate the actual current density. 
#    if inputs.flag_cathode == 1:
    i_ext_amp = -inputs.C_rate*oneC
    
    sigma_eff = inputs.sigma_cat*eps_C_0/tau**3
    
    u_Li_el = inputs.D_Li_el*eps_el_0/tau**3
    
#    D_el = inputs.D_Li_el*eps_el_0/tau**3
    D_el = inputs.D_Li_el/tau**3
    
    def get_tflag():
        return cathode.t_flag
    
    def set_tflag(value):
        cathode.t_flag = value
    
"============================================================================="        
        
class sep():
    # Set a flag to let the solver know whether to implement this class
    flag = inputs.flag_sep
    
    # Number of nodes in the y-direction
    npoints = inputs.npoints_sep
    
    # Number of variables per node
    nVars = 1 + elyte_obj.n_species
    
    H = inputs.H_elyte  # Separator thickness [m]
    
    tau = inputs.tau_sep  # Tortuosity of separator
    
    # Geometric parameters
    epsilon = inputs.epsilon_sep  # Volume fraction of separator material [-]
    epsilon_el = 1 - epsilon      # Volume fraction of electrolyte [-]
    dyInv = npoints/H             # Inverse of y-direction discretization [1/m]
    dy = H/npoints
    
    # Mobility of electrolyte species
    u_Li_el = inputs.D_Li_el*epsilon_el/ct.gas_constant/inputs.T/tau**3
    
    ptr = {}
    ptr['rho_k_el'] = np.arange(0, elyte_obj.n_species)
    ptr['phi'] = elyte_obj.n_species
    
    ptr_vec = {}
    ptr_vec['rho_k_el'] = cathode.nSV + ptr['rho_k_el']
    ptr_vec['phi'] = cathode.nSV + ptr['phi']
    
    for i in np.arange(1, npoints):
        ptr_vec['rho_k_el'] = np.append(ptr_vec['rho_k_el'], 
                                      cathode.nSV + ptr['rho_k_el'] + i*nVars)
        ptr_vec['phi'] = np.append(ptr_vec['phi'], 
                                   cathode.nSV + ptr['phi'] + i*nVars)
        
    # Set the length of the solution vector for the separator
    nSV = npoints*nVars
    
    D_el = inputs.D_Li_el*epsilon_el/tau**3
    
    offsets = np.arange(int(cathode.nSV), int(cathode.nSV) + int(nSV), int(nVars))
    
"============================================================================="

class anode():
    flag = inputs.flag_anode
    
    npoints = inputs.npoints_anode
#    
#    nshells = inputs.nshells_anode
    
    nVars = 2 + elyte_obj.n_species
    
    # Pointers
    ptr = {}
    ptr['iFar'] = elyte_obj.species_index(inputs.Li_species_elyte)
    
    ptr['rho_k_el'] = np.arange(0, elyte_obj.n_species)
    ptr['phi_dl'] = ptr['rho_k_el'][-1] + 1
    ptr['phi_ed'] = ptr['rho_k_el'][-1] + 2
    
    ptr_vec = {}
    ptr_vec['rho_k_el'] = cathode.nSV + sep.nSV + ptr['rho_k_el']
    
    for i in np.arange(1, npoints):
        ptr_vec['rho_k_el'] = np.append(ptr_vec['rho_k_el'],
                                       cathode.nSV + sep.nSV + ptr['rho_k_el'] + i*nVars)
    
    # Set length of solution vector for anode
    nSV = npoints*nVars
    offsets = np.arange(int(cathode.nSV + sep.nSV), 
                        int(cathode.nSV + sep.nSV) + int(nSV), int(nVars))
    
    # Geometric parameters
    eps_el = 1 - inputs.epsilon_an
    tau = inputs.tau_an
    r_p = inputs.r_p_cat
    dyInv = npoints/inputs.H_an
    dy = inputs.H_an/npoints
    H = inputs.H_an
    
    C_dl = inputs.C_dl_an
    A_Li = 1e3
    sigma_eff = inputs.sigma_an*inputs.epsilon_an/tau**3
    
    u_Li_el = inputs.D_Li_el*eps_el/tau**3
    
    D_el = inputs.D_Li_el*eps_el/tau**3
    
"============================================================================="

class sol_init():
    
    # Initialize solution vector 
    SV_0 = np.zeros([anode.nSV + sep.nSV + cathode.nSV])
    
    # Set up algebraic variable vector
    algvar = np.zeros_like(SV_0)
     
    # Cathode
    offsets = cathode.offsets
    ptr = cathode.ptr
    for j in np.arange(0, cathode.npoints):
        
        SV_0[offsets[j] + ptr['eps_S8']] = cathode.eps_S_0
        algvar[offsets[j] + ptr['eps_S8']] = 1
        
        SV_0[offsets[j] + ptr['eps_Li2S']] = cathode.eps_L_0
        algvar[offsets[j] + ptr['eps_Li2S']] = 1
        
        SV_0[offsets[j] + ptr['rho_k_el']] = inputs.C_k_el_0
        algvar[offsets[j] + ptr['rho_k_el']] = 1
        
        SV_0[offsets[j]+ptr['phi_dl']] = inputs.Cell_voltage - inputs.Phi_el_init
        algvar[offsets[j] + ptr['phi_dl']] = 1
                                           
        SV_0[offsets[j]+ptr['phi_ed']] = inputs.Cell_voltage
#        algvar[offsets[j] + ptr['phi_ed']] = 1
        
        SV_0[offsets[j]+ptr['np_S8']]=inputs.np_S8_init
        algvar[offsets[j] + ptr['np_S8']] = 1
        
        SV_0[offsets[j]+ptr['np_Li2S']] = inputs.np_Li2S_init
        algvar[offsets[j] + ptr['np_Li2S']] = 1
     
    # Separator
    offsets = sep.offsets
    ptr = sep.ptr
    for j in np.arange(0, sep.npoints):
        
        SV_0[offsets[j] + ptr['rho_k_el']] = inputs.C_k_el_0
        algvar[offsets[j] + ptr['rho_k_el']] = 1
        
        SV_0[offsets[j] + ptr['phi']] = inputs.Phi_el_init
     
    # Anode
    offsets = anode.offsets
    ptr = anode.ptr
    for j in np.arange(0, anode.npoints):
        SV_0[offsets[j] + ptr['rho_k_el']] = inputs.C_k_el_0
        algvar[offsets[j] + ptr['rho_k_el']] = 1
        
        SV_0[offsets[j] + ptr['phi_dl']] = inputs.Phi_an_init - inputs.Phi_el_init
        algvar[offsets[j] + ptr['phi_dl']] = 1
        
        SV_0[offsets[j] + ptr['phi_ed']] = inputs.Phi_an_init
        
    
                                           
"============================================================================="

print("Initialization check")