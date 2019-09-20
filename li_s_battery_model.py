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

def main():
    
    res_class = eval(inputs.test_type)
    
    plt.close('all')
    t_count = time.time()
        
    SV_0 = sol_init.SV_0
    SV_dot_0 = np.zeros_like(SV_0)
    t_0 = 0.
    t_f = 3600./inputs.C_rate
    algvar = sol_init.algvar
    
    atol = 1e-6; rtol = 1e-6; sim_output = 50
    
    rate_tag = str(inputs.C_rate)+"C"
    
    fig, axes = plt.subplots(sharey="row", figsize=(18,9), nrows=3, ncols = 2)
    plt.subplots_adjust(wspace = 0.15, hspace = 0.4)
    fig.text(0.15, 0.8, rate_tag, fontsize=20, bbox=dict(facecolor='white', alpha = 0.5))
    
    # Set up user function to build figures based on inputs
    
    "----------Equilibration----------"
    
    print('\nEquilibrating...')
    
    # Set external current to 0 for equilibration
#    cat.set_i_ext(0)
    
    # Create problem object
#    bat_eq = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
#    bat_eq.external_event_detection = True
#    bat_eq.algvar = algvar
    
    # Create simulation object
#    sim_eq = IDA(bat_eq)
#    sim_eq.atol = atol
#    sim_eq.rtol = rtol
#    sim_eq.verbosity = sim_output
#    sim_eq.make_consistent('IDA_YA_YDP_INIT')
#    
#    t_eq, SV_eq, SV_dot_eq = sim_eq.simulate(t_f)
    
    # Put solution into pandas dataframe with labeled columns
#    SV_eq_df = label_columns(t_eq, SV_eq, an.npoints, sep.npoints, cat.npoints)
    
    # Obtain tag strings for dataframe columns
#    tags = tag_strings(SV_eq_df)
    
#    plot_sim(tags, SV_eq_df, 'Equilibrating', 0, fig, axes)
#    print(SV_eq_df[tags['rho_el'][4:10]].iloc[-1])
#    
#    print('Done equilibrating\n')
    
    "------------Discharging-------------"
    
    print('Discharging...')
    
    # New initial conditions from previous simulation
#    SV_0 = SV_eq[-1, :]
#    SV_dot_0 = SV_dot_eq[-1, :]
    
    # Set external current
    cat.set_i_ext(cat.i_ext_amp)
    
    # Update problem instance initial conditions
    bat_ch = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
    bat_ch.external_event_detection = True
    bat_ch.algvar = algvar
    
    # Re-initialize simulation object
    sim_ch = IDA(bat_ch)
    sim_ch.atol = atol
    sim_ch.rtol = rtol
    sim_ch.verbosity = sim_output
    sim_ch.make_consistent('IDA_YA_YDP_INIT')
    
    t_ch, SV_ch, SV_dot_ch = sim_ch.simulate(t_f)
    
#    if hasattr(cathode, 'get_tflag'):
#        t_flag_ch = cathode.get_tflag
        
    SV_ch_df = label_columns(t_ch, SV_ch, an.npoints, sep.npoints, cat.npoints)
    
    # Obtain tag strings for dataframe columns
    tags = tag_strings(SV_ch_df)
    
    plot_sim(tags, SV_ch_df, 'Discharging', 1-1, fig, axes)
    
    print('Done Discharging\n')
    
    "--------Re-equilibration---------"
    
    if inputs.flag_req == 1:
        
        print('Re-equilibrating...')
        
        # New initial conditions from previous simulation
        SV_0 = SV_ch[-1, :]
        SV_dot_0 = SV_dot_ch[-1, :]
        
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
        
        plot_sim(tags, SV_req_df, 'Re-Equilibrating', 2-1, fig, axes)
    
        print('Done re-equilibrating\n')
    else:
        SV_req = SV_ch
        SV_dot_req = SV_dot_ch
        
    "-----------Charging-----------"
    
#    print('Charging...')
#    
#    SV_0 = SV_req[-1, :]
#    SV_dot_0 = SV_dot_req[-1, :]
#    
#    cat.set_i_ext(-cat.i_ext_amp)
#    
#    # Update problem instance initial conditions
#    bat_dch = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
#    bat_dch.external_event_detection = True
#    bat_dch.algvar = algvar
#    
#    # Re-initialize simulation object
#    sim_dch = IDA(bat_dch)
#    sim_dch.atol = atol
#    sim_dch.rtol = rtol
#    sim_dch.verbosity = sim_output
#    sim_dch.make_consistent('IDA_YA_YDP_INIT')
#    
#    t_dch, SV_dch, SV_dot_dch = sim_dch.simulate(t_f)
#    
##    if hasattr(cathode, 'get_tflag'):
##        t_flag_dch = cathode.get_tflag
#        
#    SV_dch_df = label_columns(t_dch, SV_dch, an.npoints, sep.npoints, cat.npoints)
#    
#    plot_sim(tags, SV_dch_df, 'Charging', 3, fig, axes)
#    
#    print('Done Charging\n')
    
    t_elapsed = time.time() - t_count
    print('t_cpu=', t_elapsed, '\n')
    
    return SV_eq_df, SV_ch_df, SV_req_df #, SV_dch_df
    
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
from li_s_battery_init import conductor_obj as conductor
from li_s_battery_init import elyte_obj as elyte
from math import pi, exp

class cc_cycling(Implicit_Problem):
    def res_fun(t, SV, SV_dot):
        
        res = np.zeros_like(SV)
        ptr = cat.ptr; F = ct.faraday; R = ct.gas_constant
        """Cathode CC boundary"""
        j = 0; offset = cat.offsets[int(j)]
        i_ext = cat.get_i_ext()
        s1 = {}
        s2 = {}
        
        # Set electronic current and ionic current boundary conditions
        i_el_p = i_ext
        i_io_p = 0
        N_io_p = 0        
        
        for j in np.arange(0, cat.npoints):
            
            i_el_m = i_el_p
            i_io_m = i_io_p
            N_io_m = N_io_p
            
            i_el_p = 0
            i_io_p = i_ext
            N_io_p = np.zeros_like(SV[offset + ptr['rho_k_el']])
            N_io_p[2] = i_ext/F
            
            np_S = SV[offset + ptr['np_S8']]
            np_L = SV[offset + ptr['np_Li2S']]
            eps_el = 1 - cat.eps_C_0 - SV[offset + ptr['eps_S8']] - SV[offset + ptr['eps_Li2S']]
            
            # Calculate new particle radii based on new volume fractions
            A_S = 3*(3*SV[offset + ptr['eps_S8']]*cat.V_0/2/pi/np_S)**(-1/3)
            A_L = 3*(3*SV[offset + ptr['eps_Li2S']]*cat.V_0/2/pi/np_L)**(-1/3)
#            A_L = cat.A_L_0*(SV[offset + ptr['eps_Li2S']]/cat.eps_L_0)**1.5
            
            r_S = 3/A_S
            r_L = 3/A_L
            
            A_C = 3*(2*cat.r_C**2 - np_S*r_S**2 - np_L*r_L**2)/2/cat.r_C**3
            
#            u_Li_el = inputs.D_Li_el*eps_el/cat.tau**3
#            D_el = inputs.D_Li_el*eps_el/cat.tau**3
                                    
            # Set states
            phi_el = SV[offset + ptr['phi_el']]
            phi_dl = SV[offset + ptr['phi_dl']]
            
            carbon.electric_potential = phi_dl + phi_el
            elyte.electric_potential = phi_el
            conductor.electric_potential = phi_dl + phi_el
            
#            rho = sum(SV[offset + ptr['rho_k_el']])
            X = SV[offset + ptr['rho_k_el']]/sum(SV[offset + ptr['rho_k_el']])

            elyte.X = X
            
            sdot_C = C_el_s.net_production_rates
            sdot_S = S_el_s.net_production_rates
            sdot_L = L_el_s.net_production_rates 
            
            # Calculate respective changes in species for each interface. This
            #   is done separately due to some species being produced/consumed
            #   at two separate interfaces - i.e. S^2- is produced at the C-el
            #   interface and consumed at the Li2S-el interface which will have
            #   different areas
            R_C = sdot_C[1:-2]*A_C
            R_S = sdot_S[1:-2]*A_S
            R_L = sdot_L[1:-2]*A_L
            
            # Net rate of formation
            R_net = R_C + R_S + R_L
                        
            i_Far = sdot_C[-2]*F*A_C/cat.dyInv
#            print(S_el_s.forward_rate_constants, '\n')
#            print(C_el_s.forward_rate_constants, '\n')
#            print(L_el_s.delta_gibbs, L_el_s.forward_rate_constants, '\n')
#            print(sdot_S[0], sdot_L[0], C_el_s.delta_gibbs, L_el_s.delta_gibbs, i_ext, '\n')
            
            """Calculate change in Sulfur"""
            res[offset + ptr['eps_S8']] = (SV_dot[offset + ptr['eps_S8']] - sulfur.volume_mole*sdot_S[0]*A_S)
       
            """Calculate change in Li2S"""
            res[offset + ptr['eps_Li2S']] = (SV_dot[offset + ptr['eps_Li2S']] - Li2S.volume_mole*abs(sdot_L[0])*A_L)
            
            """Calculate change in electrolyte"""
            res[offset + ptr['rho_k_el']] = (SV_dot[offset + ptr['rho_k_el']]
            - (R_net + (N_io_m - N_io_p)*cat.dyInv)/eps_el
            + SV[offset + ptr['rho_k_el']]*(- SV_dot[offset + ptr['eps_S8']] 
                                            - SV_dot[offset + ptr['eps_Li2S']])/eps_el)
            
            """Calculate change in delta-phi double layer"""
            res[offset + ptr['phi_dl']] = SV_dot[offset + ptr['phi_dl']] - (-i_Far + i_el_m - i_el_p)*cat.dyInv/cat.C_dl/A_C
            
            """Algebraic expression for charge neutrality in all phases"""
            res[offset + ptr['phi_el']] = SV[offset + ptr['phi_el']]  #i_el_m - i_el_p + i_io_m - i_io_p
            
            """Calculate change in S8 nucleation sites"""
            res[offset + ptr['np_S8']] = SV_dot[offset + ptr['np_S8']]
            
            """Calculate change in Li2S nucleation sites"""
            res[offset + ptr['np_Li2S']] = SV_dot[offset + ptr['np_Li2S']]
            
#        print(res, '\n')
#        print(t, i_ext)
        
        """==================Separator boundary conditions=================="""
            
#        i_el_m = i_el_p
#        i_io_m = i_io_p
#        N_io_m = N_io_p
            
        # Set ionic current and flux at the separator boundary
#        i_el_p = 0
#        i_io_p = i_ext
#        N_io_p = i_ext/F
        
        return res  
      
    "========================================================================="
    
    def set_state():
        return
    
    "========================================================================="
    
    def state_events(self, t, y, yd, sw):
        
        event1 = np.zeros([cat.npoints*int(y[cat.ptr_vec['np_S8']])])
        event2 = np.zeros([cat.npoints*int(y[cat.ptr_vec['np_S8']])])
        event1 = 1 - y[cat.ptr_vec['eps_S8']]
        event2 = y[cat.ptr_vec['eps_S8']]
        
        event3 = np.zeros([cat.npoints*int(y[cat.ptr_vec['np_Li2S']])])
        event4 = np.zeros([cat.npoints*int(y[cat.ptr_vec['np_Li2S']])])
        event3 = 1 - y[cat.ptr_vec['eps_Li2S']]
        event4 = y[cat.ptr_vec['eps_Li2S']]
        
        events = np.concatenate((event1, event2, event3, event4))

        return events
    
    "========================================================================="
    
    def handle_event(self, solver, event_info):
        
        state_info = event_info[0]
        
        if any(state_info):
            print('shoulda stopped')
            raise TerminateSimulation
#        while True:
#            self.event_switch(solver, event_info)
#            self.init_mode(solver)
            
#            if not True in event_info:
#                break
    
    "========================================================================="
    
    def event_switch(self, solver, event_info):
        if not all(event_info):
            solver.sw = [not solver.sw]
            
        return solver.sw
    
    "========================================================================="
    
    def init_mode(self, solver):
        cat.set_tflag(solver.t)
        solver.make_consistent('IDA_YA_YDP_INIT')
        
        if battery.get_i_ext() != 0:
            battery.set_i_ext(0)
    
    "========================================================================="
    
#class is_cycling(Implicit_Problem):
#    def res_fun():
#        print('is_cycling')
    
    
if __name__ == "__main__":
    SV_eq_df, SV_ch_df, SV_req_df = main()
#    SV_eq, SV_ch, SV_req, SV_dch = main()

