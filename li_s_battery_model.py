# -*- coding: utf-8 -*-
"""
Created on Thu May  9 12:30:28 2019

@author: dkorff

This module contains the set up and running of the simulation for the Li-S
model.
"""

import numpy as np
import time
import importlib
import cantera as ct
from matplotlib import pyplot as plt

from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem
from assimulo.exception import TerminateSimulation

from li_s_battery_inputs import inputs

#from li_s_battery_init import battery
from li_s_battery_init import anode as an
from li_s_battery_init import sep
from li_s_battery_init import cathode as cat
from li_s_battery_init import sol_init

from li_s_battery_post import label_columns
from li_s_battery_post import tag_strings
from li_s_battery_post import plot_sim
from li_s_battery_post import plot_meanPS

def main():
    
    res_class = eval(inputs.test_type)
    
#    plt.close('all')
    t_count = time.time()
        
    SV_0 = sol_init.SV_0
    SV_dot_0 = np.zeros_like(SV_0)
    t_0 = 0.
    t_f = 3600./inputs.C_rate
    algvar = sol_init.algvar
    atol = np.ones_like(SV_0)*1e-5
    atol[cat.ptr_vec['eps_S8']] = 1e-30
    atol[cat.ptr_vec['eps_Li2S']] = 1e-25
    atol[cat.ptr_vec['rho_k_el']] = 1e-30
#    atol = 1e-30; 
    rtol = 1e-6; sim_output = 50
    
    rate_tag = str(inputs.C_rate)+"C"
    
    fig, axes = plt.subplots(sharey="row", figsize=(9,12), nrows=3, ncols = 2+inputs.flag_req)
    plt.subplots_adjust(wspace = 0.15, hspace = 0.4)
    fig.text(0.35, 0.85, rate_tag, fontsize=20, bbox=dict(facecolor='white', alpha = 0.5))
    
#    fig2, axes2 = plt.subplots(sharey="row", figsize=(9,12), nrows=3, ncols = 1)
#    plt.subplots_adjust(wspace = 0.15, hspace = 0.4)
#    fig.text(0.15, 0.8, rate_tag, fontsize=20, bbox=dict(facecolor='white', alpha = 0.5))
    
    # Set up user function to build figures based on inputs
    
    "----------Equilibration----------"
    
    print('\nEquilibrating...')
    
    # Set external current to 0 for equilibration
    cat.set_i_ext(0)
    
    # Create problem object
    bat_eq = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
    bat_eq.external_event_detection = True
    bat_eq.algvar = algvar
    
    # Create simulation object
    sim_eq = IDA(bat_eq)
    sim_eq.atol = atol
    sim_eq.rtol = rtol
    sim_eq.verbosity = sim_output
    sim_eq.make_consistent('IDA_YA_YDP_INIT')
    
    t_eq, SV_eq, SV_dot_eq = sim_eq.simulate(t_f)
    
    # Put solution into pandas dataframe with labeled columns
    SV_eq_df = label_columns(t_eq, SV_eq, an.npoints, sep.npoints, cat.npoints)
#    SV_eq_df = []
    
    # Obtain tag strings for dataframe columns
    tags = tag_strings(SV_eq_df)
    
#    plot_sim(tags, SV_eq_df, 'Equilibrating', 0, fig, axes)
#    print(SV_eq_df[tags['rho_el'][4:10]].iloc[-1])
    
    print('Done equilibrating\n')
    
    "------------Discharging-------------"
    
    print('Discharging...')
    
    # New initial conditions from previous simulation
    SV_0 = SV_0  #SV_eq[-1, :]
    SV_dot_0 = SV_dot_0  #SV_dot_eq[-1, :]
    
    # Set external current
    cat.set_i_ext(cat.i_ext_amp)
    
    # Update problem instance initial conditions
    bat_dch = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
    bat_dch.external_event_detection = True
    bat_dch.algvar = algvar
        
    # Re-initialize simulation object
    sim_dch = IDA(bat_dch)
    sim_dch.atol = atol
    sim_dch.rtol = rtol
#    sim_dch.maxh = 5
    sim_dch.verbosity = sim_output
    sim_dch.make_consistent('IDA_YA_YDP_INIT')
    
    t_dch, SV_dch, SV_dot_dch = sim_dch.simulate(t_f)
    
#    if hasattr(cathode, 'get_tflag'):
#        t_flag_ch = cathode.get_tflag
        
    SV_dch_df = label_columns(t_dch, SV_dch, an.npoints, sep.npoints, cat.npoints)
#    SV_dch_df = []
    # Obtain tag strings for dataframe columns
    tags = tag_strings(SV_dch_df)
    
    plot_sim(tags, SV_dch_df, 'Discharging', 0, fig, axes)
    
    plot_meanPS(SV_dch_df, tags, 'Discharging')
    
    print('Done Discharging\n')
    
    "--------Re-equilibration---------"
    
    if inputs.flag_req == 1:
        
        print('Re-equilibrating...')
        
        # New initial conditions from previous simulation
        SV_0 = SV_dch[-1, :]
        SV_dot_0 = SV_dot_dch[-1, :]
        
        # Set external current
        cat.set_i_ext(0)
        
        # Update problem instance initial conditions
        bat_req = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
        bat_req.external_event_detection = True
        bat_req.algvar = algvar
        
        # Re-initialize simulation object
        sim_req = IDA(bat_req)
        sim_req.atol = atol
        sim_req.rtol = rtol
        sim_req.verbosity = sim_output
        sim_req.make_consistent('IDA_YA_YDP_INIT')
        
        t_req, SV_req, SV_dot_req = sim_req.simulate(t_f)
        
        SV_req_df = label_columns(t_req, SV_req, an.npoints, sep.npoints, cat.npoints)
        
        plot_sim(tags, SV_req_df, 'Re-Equilibrating', 1, fig, axes)
    
        print('Done re-equilibrating\n')
    else:
        SV_req = SV_dch
        SV_dot_req = SV_dot_dch
        
    "-----------Charging-----------"
    
    print('Charging...')
    
    SV_0 = SV_req[-1, :]  #SV_dch[-1, :]
    SV_dot_0 = SV_dot_req[-1, :]  #SV_dot_dch[-1, :]
    
    SV_0[0] = cat.eps_cutoff  # cat.eps_cutoff*1e2
#    SV_dot_0[0] = 0
    
    cat.set_i_ext(-cat.i_ext_amp)
    
    # Update problem instance initial conditions
    bat_ch = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
    bat_ch.external_event_detection = True
    bat_ch.algvar = algvar
    
    # Re-initialize simulation object
    sim_ch = IDA(bat_ch)
    sim_ch.atol = atol
    sim_ch.rtol = rtol
#    sim_ch.maxh = 0.1
    sim_ch.verbosity = sim_output
    sim_ch.make_consistent('IDA_YA_YDP_INIT')
    
    t_ch, SV_ch, SV_dot_ch = sim_ch.simulate(t_f)
    
#    if hasattr(cathode, 'get_tflag'):
#        t_flag_dch = cathode.get_tflag
        
    SV_ch_df = label_columns(t_ch, SV_ch, an.npoints, sep.npoints, cat.npoints)
    
    plot_sim(tags, SV_ch_df, 'Charging', 1+inputs.flag_req, fig, axes)
    
    plot_meanPS(SV_ch_df, tags, 'Charging')
    
    print('Max S_8(e) concentration = ', max(SV_ch[:, 6]))
    
    print('Done Charging\n')
    
    t_elapsed = time.time() - t_count
    print('t_cpu=', t_elapsed, '\n')
    
    return SV_eq_df, SV_dch_df, SV_ch_df, tags # SV_ch_df, tags #SV_eq_df, SV_req_df #, SV_dch_df
    
"=============================================================================" 
"===========RESIDUAL CLASSES AND HELPER FUNCTIONS BEYOND THIS POINT==========="
"============================================================================="

#from li_s_battery_init import anode_obj as anode
#from li_s_battery_init import anode_s_obj as anode_s
from li_s_battery_init import sulfur_obj as sulfur
from li_s_battery_init import Li2S_obj as Li2S
from li_s_battery_init import carbon_obj as carbon
from li_s_battery_init import sulfur_el_s as S_el_s
from li_s_battery_init import Li2S_el_s as L_el_s
from li_s_battery_init import carbon_el_s as C_el_s
from li_s_battery_init import Li2S_tpb
from li_s_battery_init import lithium_obj as lithium
from li_s_battery_init import lithium_el_s as lithium_s
from li_s_battery_init import conductor_obj as conductor
from li_s_battery_init import elyte_obj as elyte
from li_s_battery_functions import set_state
from li_s_battery_functions import set_state_sep
from li_s_battery_functions import dst
from math import pi, exp, tanh

class cc_cycling(Implicit_Problem):    
    def res_fun(t, SV, SV_dot):
        
        res = np.zeros_like(SV)
        ptr = cat.ptr; F = ct.faraday; R = ct.gas_constant; T = inputs.T
        
        """=============================CATHODE============================="""
        """CC BOUNDARY"""
        j = 0; offset = cat.offsets[int(j)]
        i_ext = cat.get_i_ext()
        s2 = set_state(SV, offset, cat.ptr)

        # Set electronic current and ionic current boundary conditions
        i_el_p = i_ext
        i_io_p = 0
        N_io_p = 0  
        
        """=============================CATHODE============================="""
        """INTERIOR NODES"""
        for j in np.arange(1, cat.npoints):
            
            # Set previous outlet fluxes to new inlet fluxes
            i_el_m = i_el_p
            i_io_m = i_io_p
            N_io_m = N_io_p
            s1 = dict(s2)
            
            # Update offset to NEXT node
            offset = cat.offsets[int(j)]
            
            s2 = set_state(SV, offset, cat.ptr)
            
            # Set variables to CURRENT NODE value
            offset = cat.offsets[int(j-1)]
            np_S = SV[offset + ptr['np_S8']]
            np_L = SV[offset + ptr['np_Li2S']]
            eps_S8 = max(SV[offset + ptr['eps_S8']], cat.eps_cutoff)
            eps_Li2S = max(SV[offset + ptr['eps_Li2S']], cat.eps_cutoff)
            eps_el = 1 - cat.eps_C_0 - eps_S8 - eps_Li2S  
            
            # Calculate new particle radii based on new volume fractions
            A_S = 3*eps_S8**(2/3)/(3/2/pi/np_S)**(1/3)
            A_L = 3*eps_Li2S/(3*eps_Li2S/2/pi/np_L)**(1/3)
            
            r_S = 3*eps_S8/A_S
            r_L = 3*eps_Li2S/A_L
            
            A_C = inputs.A_C_0 - (pi*np_S*r_S**2) - (pi*np_L*r_L**2)

            # Set states for THIS node            
            carbon.electric_potential = s1['phi_ed']
            elyte.electric_potential = s1['phi_el'] 
            conductor.electric_potential = s1['phi_ed']
            elyte.X = s1['X_k']
            
            D_el = cat.D_el*eps_el**(1.5)

            # Current node plus face boundary fluxes
            i_el_p = cat.sigma_eff*(s1['phi_ed'] - s2['phi_ed'])*cat.dyInv
            N_io_p, i_io_p = dst(s1, s2, D_el, cat.dy, cat.dy)
            
            sdot_C = C_el_s.get_net_production_rates(elyte)
            sdot_L = L_el_s.get_net_production_rates(elyte) 
            
            # Calculate respective changes in species for each interface. This
            #   is done separately due to some species being produced/consumed
            #   at two separate interfaces - i.e. S^2- is produced at the C-el
            #   interface and consumed at the Li2S-el interface which will have
            #   different areas
            R_C = sdot_C*A_C
            R_L = sdot_L*A_L
            if SV[offset + ptr['eps_S8']] < cat.eps_cutoff:
                sdot_S = 0*S_el_s.get_net_production_rates(elyte)
                R_S = sdot_S*A_S
            else:
                sdot_S = S_el_s.get_net_production_rates(elyte)
                R_S = sdot_S*A_S
                
            i_Far = C_el_s.get_net_production_rates(conductor)*F*A_C/cat.dyInv
            
            # Net rate of formation
            R_net = R_C + R_S + R_L
            R_net[cat.ptr['iFar']] += (-i_Far + i_el_m - i_el_p)/cat.dy/F
            
            sdot_S8 = S_el_s.get_net_production_rates(sulfur)
            sdot_Li2S = L_el_s.get_net_production_rates(Li2S)  
            
            """Calculate change in Sulfur"""                
            res[offset + ptr['eps_S8']] = (SV_dot[offset + ptr['eps_S8']] 
                                        - sulfur.volume_mole*sdot_S8*A_S)
       
            """Calculate change in Li2S"""
            res[offset + ptr['eps_Li2S']] = (SV_dot[offset + ptr['eps_Li2S']] 
                                          - Li2S.volume_mole*sdot_Li2S*A_L)
            
            """Calculate change in electrolyte"""
            res[offset + ptr['rho_k_el']] = (SV_dot[offset + ptr['rho_k_el']] - 
            (R_net + (N_io_m - N_io_p)*cat.dyInv)/eps_el 
            + SV[offset + ptr['rho_k_el']]*(- SV_dot[offset + ptr['eps_S8']] 
                                            - SV_dot[offset + ptr['eps_Li2S']])/eps_el)
            
            """Calculate change in delta-phi double layer"""
            res[offset + ptr['phi_dl']] = (SV_dot[offset + ptr['phi_dl']] - 
            (-i_Far + i_el_m - i_el_p)*cat.dyInv/cat.C_dl/A_C)
            
            """Algebraic expression for charge neutrality in all phases"""
            res[offset + ptr['phi_ed']] = i_el_m - i_el_p + i_io_m - i_io_p
            
            """Calculate change in S8 nucleation sites"""
            res[offset + ptr['np_S8']] = SV_dot[offset + ptr['np_S8']]
            
            """Calculate change in Li2S nucleation sites"""
            res[offset + ptr['np_Li2S']] = SV_dot[offset + ptr['np_Li2S']]
            
        """=============================CATHODE============================="""
        """SEPARATOR BOUNDARY"""
        
        i_el_m = i_el_p
        i_io_m = i_io_p
        N_io_m = N_io_p
        s1 = dict(s2)
        
        # Shift forward to NEXT node, first separator node (j=0)
        j = 0; offset = sep.offsets[int(j)]
        
        s2 = set_state_sep(SV, offset, sep.ptr)
        
        # Shift back to THIS node, set THIS node outlet conditions
        j = cat.npoints-1; offset = cat.offsets[int(j)]
        
        # Set variables to CURRENT NODE value
        np_S = SV[offset + ptr['np_S8']]
        np_L = SV[offset + ptr['np_Li2S']]
        eps_S8 = max(SV[offset + ptr['eps_S8']], cat.eps_cutoff)
        eps_Li2S = max(SV[offset + ptr['eps_Li2S']], cat.eps_cutoff)
        eps_el = 1 - cat.eps_C_0 - eps_S8 - eps_Li2S  
        
        # Calculate new particle radii based on new volume fractions
        A_S = 2*pi*np_S*(3*eps_S8/2/np_S/pi)**(2/3)
        A_L = 2*pi*np_L*(3*eps_Li2S/2/np_L/pi)**(2/3)
        
        r_S = 3*eps_S8/A_S
        r_L = 3*eps_Li2S/A_L
        
        tpb_len = 3*eps_Li2S/(r_L**2)
        
        A_C = inputs.A_C_0 - (pi*np_S*r_S**2) - (pi*np_L*r_L**2)
        
        carbon.electric_potential = s1['phi_ed']
        elyte.electric_potential = s1['phi_el']
        conductor.electric_potential = s1['phi_ed']
        
        elyte.X = s1['X_k']
        
        # Set outlet boundary conditions for THIS node
        i_el_p = 0
        D_el = cat.D_el*eps_el**(1.5)
        dyInv_boundary = 1/(0.5*(cat.dy + sep.dy))
        N_io_p, i_io_p = dst(s1, s2, D_el, cat.dy, sep.dy)
        
        sdot_C = C_el_s.get_net_production_rates(elyte)
#        sdot_L = L_el_s.get_net_production_rates(elyte)
#        print(C_el_s.forward_rate_constants, '\n\n', C_el_s.reverse_rate_constants, i_ext, '\n\n\n')
#        R_tpb = tpb_len*Li2S_tpb.get_net_production_rates(elyte)
        R_C = sdot_C*A_C
#        R_L = sdot_L*A_L
        if eps_S8 < cat.eps_dropoff and i_ext > 0:  #eps_S8 < cat.eps_dropoff and i_ext > 0: #S_el_s.get_net_production_rates(sulfur) < 0 and
            mult = (1/cat.eps_dropoff)*tanh(eps_S8)
            sdot_S8 = S_el_s.get_creation_rates(sulfur) - mult*S_el_s.get_destruction_rates(sulfur)
            sdot_S = S_el_s.get_net_production_rates(elyte)
            R_S = sdot_S*A_S
        else:
            sdot_S8 = S_el_s.get_net_production_rates(sulfur)
            sdot_S = S_el_s.get_net_production_rates(elyte)
            R_S = sdot_S*A_S
                        
        if eps_Li2S < cat.eps_dropoff and i_ext < 0:
#            R_tpb = 0*tpb_len*Li2S_tpb.get_net_production_rates(elyte)
            mult = (1/cat.eps_dropoff)*tanh(eps_Li2S)
            sdot_Li2S = L_el_s.get_creation_rates(Li2S) - mult*(L_el_s.get_destruction_rates(Li2S))
            sdot_L = 0*L_el_s.get_net_production_rates(elyte)
            R_L = sdot_L*A_L
        else:
#            R_tpb = tpb_len*Li2S_tpb.get_net_production_rates(elyte)
            sdot_Li2S = L_el_s.get_net_production_rates(Li2S)
            sdot_L = L_el_s.get_net_production_rates(elyte)
            R_L = sdot_L*A_L
             
        i_C = C_el_s.get_net_production_rates(conductor)*A_C
        i_L = 0*Li2S_tpb.get_net_production_rates(conductor)*tpb_len
        i_Far = (i_C + i_L)*F/cat.dyInv
        
        # Net rate of formation
        R_net = R_C + R_S + R_L #+ R_tpb
        R_net[cat.ptr['iFar']] += (-i_Far + i_el_m - i_el_p)/cat.dy/F
        
#        sdot_S8 = S_el_s.get_net_production_rates(sulfur)
#        sdot_Li2S = L_el_s.get_net_production_rates(Li2S)  
        sdot_tpb = Li2S_tpb.get_net_production_rates(Li2S)
                        
        """Calculate change in Sulfur"""                
        res[offset + ptr['eps_S8']] = (SV_dot[offset + ptr['eps_S8']] 
                                    - sulfur.volume_mole*sdot_S8*A_S)
        
        """Calculate change in Li2S"""
        res[offset + ptr['eps_Li2S']] = (SV_dot[offset + ptr['eps_Li2S']] 
                                      - Li2S.volume_mole*sdot_Li2S*A_L) # - Li2S.volume_mole*sdot_tpb*tpb_len)
                
        """Calculate change in electrolyte"""
        res[offset + ptr['rho_k_el']] = (SV_dot[offset + ptr['rho_k_el']] - 
        (R_net + (N_io_m - N_io_p)*cat.dyInv)/eps_el
        + SV[offset + ptr['rho_k_el']]*(- SV_dot[offset + ptr['eps_S8']] 
                                        - SV_dot[offset + ptr['eps_Li2S']])/eps_el)
        
        """Calculate change in delta-phi double layer"""
        res[offset + ptr['phi_dl']] = (SV_dot[offset + ptr['phi_dl']] - 
        (-i_Far + i_el_m - i_el_p)*cat.dyInv/cat.C_dl/A_C)
        
        """Algebraic expression for charge neutrality in all phases"""
        res[offset + ptr['phi_ed']] = i_el_m - i_el_p + i_io_m - i_io_p
        
        """Calculate change in S8 nucleation sites"""
        res[offset + ptr['np_S8']] = SV_dot[offset + ptr['np_S8']]
        
        """Calculate change in Li2S nucleation sites"""
        res[offset + ptr['np_Li2S']] = SV_dot[offset + ptr['np_Li2S']]
        
        """============================SEPARATOR============================"""
        """INTERIOR NODES"""
        
        """============================SEPARATOR============================"""
        """CATHODE BOUNDARY"""
        
        i_io_m = i_io_p
        N_io_m = N_io_p
        s1 = dict(s2)
        
        # Shift forward to NEXT node
        j = 0; offset = an.offsets[int(j)]
        s2 = set_state(SV, offset, an.ptr)
                
        # Shift back to THIS node
        offset = sep.offsets[int(j)]
        
        D_el = sep.D_el 
        
        # Current node plus face boundary conditions
        N_io_p, i_io_p = dst(s1, s2, D_el, sep.dy, an.dy)
        
        res[offset + sep.ptr['rho_k_el']] = (SV_dot[offset + sep.ptr['rho_k_el']]
        - (N_io_m - N_io_p)*sep.dyInv/sep.epsilon_el)
                
        res[offset + sep.ptr['phi']] = i_io_m - i_io_p
        
        """==============================ANODE=============================="""
        """INTERIOR NODES"""
          
        i_io_m = i_io_p
        N_io_m = N_io_p
        i_el_m = 0
        s1 = dict(s2)
        
        j = 0
        offset = an.offsets[int(j)]
        
        i_el_p = i_ext
        i_io_p = 0
        N_io_p = 0
        
        elyte.X = s1['X_k']
        elyte.electric_potential = s1['phi_el']
        lithium.electric_potential = s1['phi_ed']
        conductor.electric_potential = s1['phi_ed']
        
        sdot_Li = lithium_s.get_net_production_rates(elyte)
        sdot_Far = lithium_s.get_net_production_rates(conductor)
        
        R_net = sdot_Li*an.A_Li
        i_Far = sdot_Far*an.A_Li*F*an.dy
        
        res[offset + an.ptr['rho_k_el']] = (SV_dot[offset + an.ptr['rho_k_el']]
        - (R_net + (N_io_m - N_io_p)*an.dyInv)/an.eps_el)

        res[offset + an.ptr['phi_dl']] = (SV_dot[offset + an.ptr['phi_dl']]
        - (-i_Far + i_el_m - i_el_p)*an.dyInv/an.C_dl/an.A_Li) 
        
        res[offset + an.ptr['phi_ed']] = SV[offset + an.ptr['phi_ed']]
        
        """==============================ANODE=============================="""
        """CC BOUNDARY"""
#        print(SV, '\n')
#        print(res, t, '\n\n')
#        print(SV_dot, '\n\n')
#        print(A_S, t, i_ext)
#        print(t)
#        if i_ext > 0:
#            print(res, '\n\n')
        
        return res  
    
    "========================================================================="
    
    def state_events(self, t, y, yd, sw):
        
        event1 = np.zeros([cat.npoints])
        event2 = np.zeros([cat.npoints])
        event1 = 1 - y[cat.ptr_vec['eps_S8']]
#        event2 = y[cat.ptr_vec['eps_S8']]
        
        event3 = np.zeros([cat.npoints])
        event4 = np.zeros([cat.npoints])
        event3 = 1 - y[cat.ptr_vec['eps_Li2S']]
#        event4 = y[cat.ptr_vec['eps_Li2S']]
        
        event5 = np.zeros([cat.npoints])
        event5 = 2.8 - y[cat.ptr_vec['phi_ed']]
        event6 = np.zeros([cat.npoints])
        event6 = y[cat.ptr_vec['phi_ed']] - 1.5
        
        event7 = np.zeros([cat.npoints*elyte.n_species])
#        event7 = y[cat.ptr_vec['rho_k_el']] - 1e-20
        
        
        events = np.concatenate((event1, event2, event3, event4, event5, event6,
                                 event7))

        return events
    
    "========================================================================="
    
    def handle_event(self, solver, event_info):
        
        state_info = event_info[0]
        
        if state_info[0]:
            print('Sulfur volume fraction over 1')
            raise TerminateSimulation
        if state_info[1]:
            print('WARNING: Sulfur volume fraction below 0')
#            raise TerminateSimulation
        elif state_info[2]:
            print('Li2S volume fraction over 1')
            raise TerminateSimulation
        elif state_info[3]:
            print('Li2S volume fraction below 0')
            raise TerminateSimulation
        elif state_info[4]:
            print('Cell voltage hit 2.8')
            raise TerminateSimulation
        elif state_info[5]:
            print('Cell voltage hit 1.5')
            raise TerminateSimulation
        elif any(state_info):
            print('Stop condition')
            raise TerminateSimulation
#        if any(state_info):
#            print('shoulda stopped')
#            raise TerminateSimulation
#        while True:
#            self.event_switch(solver, event_info)
#            self.init_mode(solver)
#            
#            if not True in event_info:
#                break
    
    "========================================================================="
    
#    def event_switch(self, solver, event_info):
#        if not all(event_info):
#            solver.sw = [not solver.sw]
#            
#        return solver.sw
    
    "========================================================================="
    
#    def init_mode(self, solver):
#        cat.set_tflag(solver.t)
#        solver.make_consistent('IDA_YA_YDP_INIT')
#        
#        if battery.get_i_ext() != 0:
#            battery.set_i_ext(0)
    
    "========================================================================="
    
#class is_cycling(Implicit_Problem):
#    def res_fun():
#        print('is_cycling')
    
    
if __name__ == "__main__":
#    SV_eq, SV_dch, tags = main()
    SV_eq, SV_dch, SV_ch, tags = main()
#    SV_eq_df, SV_ch_df, SV_req_df = main()
#    SV_eq, SV_ch, SV_req, SV_dch = main()

