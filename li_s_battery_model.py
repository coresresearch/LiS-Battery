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
    atol = np.ones_like(SV_0)*1e-6
    atol[cat.ptr_vec['eps_S8']] = 1e-25
    atol[cat.ptr_vec['eps_Li2S']] = 1e-20
    atol[cat.ptr_vec['rho_k_el']] = 1e-25
#    atol = 1e-30; 
    rtol = 1e-5; sim_output = 50
    
    rate_tag = str(inputs.C_rate)+"C"
    
    fig, axes = plt.subplots(sharey="row", figsize=(9,12), nrows=3, ncols = 1)
    plt.subplots_adjust(wspace = 0.15, hspace = 0.4)
    fig.text(0.15, 0.8, rate_tag, fontsize=20, bbox=dict(facecolor='white', alpha = 0.5))
    
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
    SV_0 = SV_eq[-1, :]
    SV_dot_0 = SV_dot_eq[-1, :]
    
    # Set external current
    cat.set_i_ext(cat.i_ext_amp)
    cc_cycling.set_Q_dl(0)
    
    # Update problem instance initial conditions
    bat_dch = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
    bat_dch.external_event_detection = True
    bat_dch.algvar = algvar
    
    # Re-initialize simulation object
    sim_dch = IDA(bat_dch)
    sim_dch.atol = atol
    sim_dch.rtol = rtol
    sim_dch.verbosity = sim_output
    sim_dch.make_consistent('IDA_YA_YDP_INIT')
    
    t_dch, SV_dch, SV_dot_dch = sim_dch.simulate(t_f)
    
#    if hasattr(cathode, 'get_tflag'):
#        t_flag_ch = cathode.get_tflag
        
    SV_dch_df = label_columns(t_dch, SV_dch, an.npoints, sep.npoints, cat.npoints)
#    SV_dch_df = []
    # Obtain tag strings for dataframe columns
#    tags = tag_strings(SV_ch_df)
    
    plot_sim(tags, SV_dch_df, 'Discharging', 1, fig, axes)
    
    plot_meanPS(SV_dch_df, tags)
    
    print('Done Discharging\n')
    
    "--------Re-equilibration---------"
    
#    if inputs.flag_req == 1:
#        
#        print('Re-equilibrating...')
#        
#        # New initial conditions from previous simulation
#        SV_0 = SV_ch[-1, :]
#        SV_dot_0 = SV_dot_ch[-1, :]
#        
#        # Set external current
#        cat.set_i_ext(0)
#        
#        # Update problem instance initial conditions
#        bat_req = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
#        bat_req.external_event_detection = True
#        bat_req.algvar = algvar
#        
#        # Re-initialize simulation object
#        sim_req = IDA(bat_req)
#        sim_req.atol = atol
#        sim_req.rtol = rtol
#        sim_req.verbosity = sim_output
#        sim_req.make_consistent('IDA_YA_YDP_INIT')
#        
#        t_req, SV_req, SV_dot_req = sim_req.simulate(t_f)
#        
#        SV_req_df = label_columns(t_req, SV_req, an.npoints, sep.npoints, cat.npoints)
#        
##        plot_sim(tags, SV_req_df, 'Re-Equilibrating', 2-1, fig, axes)
#    
#        print('Done re-equilibrating\n')
#    else:
#        SV_req = SV_ch
#        SV_dot_req = SV_dot_ch
        
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
    
    return SV_eq_df, SV_dch_df, tags #SV_eq_df, SV_req_df #, SV_dch_df
    
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
from li_s_battery_init import lithium_obj as lithium
from li_s_battery_init import lithium_el_s as lithium_s
from li_s_battery_init import conductor_obj as conductor
from li_s_battery_init import elyte_obj as elyte
from li_s_battery_functions import set_state
from li_s_battery_functions import set_state_sep
from li_s_battery_functions import dst
from math import pi

class cc_cycling(Implicit_Problem):
    
    
    def get_Q_dl():
        return cc_cycling.Q_dl
    
    def set_Q_dl(value):
        cc_cycling.Q_dl = value
    
    def res_fun(t, SV, SV_dot):
        
        res = np.zeros_like(SV)
        ptr = cat.ptr; F = ct.faraday; R = ct.gas_constant; T = inputs.T
        """Cathode CC boundary"""
        j = 0; offset = cat.offsets[int(j)]
        i_ext = cat.get_i_ext()
        s2 = {}
        
        # Set electronic current and ionic current boundary conditions
        i_el_p = i_ext
        i_io_p = 0
        N_io_p = 0  
        
        for j in np.arange(0, cat.npoints):
            
            # Set previous outlet fluxes to new inlet fluxes
            i_el_m = i_el_p
            i_io_m = i_io_p
            N_io_m = N_io_p
            s1 = dict(s2)
            
            # Update offset to current node
            offset = cat.offsets[int(j)]
            
            # Set variables to loop value
            np_S = SV[offset + ptr['np_S8']]
            np_L = SV[offset + ptr['np_Li2S']]
            eps_S8 = max(SV[offset + ptr['eps_S8']], 1e-25)
            eps_Li2S = max(SV[offset + ptr['eps_Li2S']], 1e-25)
            eps_el = 1 - cat.eps_C_0 - eps_S8 - eps_Li2S  
            
            
            # Calculate new particle radii based on new volume fractions
            A_S = 3*eps_S8/(3*eps_S8*cat.V_0/2/pi/np_S)**(1/3)
            A_L = 3*eps_Li2S/(3*eps_Li2S*cat.V_0/2/pi/np_L)**(1/3)
            
            r_S = 3*eps_S8/A_S
            r_L = 3*eps_Li2S/A_L
            
            A_C = inputs.A_C_0 - (pi*np_S*r_S**2)/cat.V_0 - (pi*np_L*r_L**2)/cat.V_0

            # Set states for THIS node
            s1 = set_state(SV, offset, cat.ptr)
            
            phi_ed = s1['phi_ed']
            phi_el = s1['phi_el']
            
            carbon.electric_potential = phi_ed
            elyte.electric_potential = phi_el 
            conductor.electric_potential = phi_ed
            
            elyte.X = s1['X_k']
            
            # Shift forward to NEXT node
            offset = sep.offsets[int(j)]
            s2 = set_state_sep(SV, offset, sep.ptr)
            
            dyInv_boundary = 1/(0.5*(cat.dy + sep.dy))
            
            # Shift back to THIS node
            offset = cat.offsets[int(j)]
            
            D_el = cat.D_el*eps_el

            # Current node plus face boundary fluxes
            i_el_p = 0
            N_io_p, i_io_p = dst(s1, s2, D_el, dyInv_boundary)
            
            sdot_C = C_el_s.get_net_production_rates(elyte)
            sdot_L = L_el_s.get_net_production_rates(elyte) 
            
            # Calculate respective changes in species for each interface. This
            #   is done separately due to some species being produced/consumed
            #   at two separate interfaces - i.e. S^2- is produced at the C-el
            #   interface and consumed at the Li2S-el interface which will have
            #   different areas
            R_C = sdot_C*A_C
            R_L = sdot_L*A_L
            if SV[offset + ptr['eps_S8']] < 0:
                sdot_S = S_el_s.get_net_production_rates(elyte)
                R_S = 0*sdot_S*A_S
            else:
                sdot_S = S_el_s.get_net_production_rates(elyte)
                R_S = sdot_S*A_S
            
            i_Far = C_el_s.get_net_production_rates(conductor)*F*A_C/cat.dyInv
            
            # Net rate of formation
            R_net = R_C + R_S + R_L
            R_net[cat.ptr['iFar']] += (-i_Far + i_el_m - i_el_p)*A_C/F
            
            sdot_S8 = S_el_s.get_net_production_rates(sulfur)
            sdot_Li2S = L_el_s.get_net_production_rates(Li2S)
            
            """Calculate change in Sulfur"""                
            res[offset + ptr['eps_S8']] = (SV_dot[offset + ptr['eps_S8']] - sulfur.volume_mole*sdot_S8*A_S)
       
            """Calculate change in Li2S"""
            res[offset + ptr['eps_Li2S']] = (SV_dot[offset + ptr['eps_Li2S']] - Li2S.volume_mole*sdot_Li2S*A_L)
            
            """Calculate change in electrolyte"""
            res[offset + ptr['rho_k_el']] = (SV_dot[offset + ptr['rho_k_el']]
            - (R_net + (N_io_m - N_io_p)*cat.dyInv)/eps_el
            + SV[offset + ptr['rho_k_el']]*(- SV_dot[offset + ptr['eps_S8']] 
                                            - SV_dot[offset + ptr['eps_Li2S']])/eps_el)
            
            """Calculate change in delta-phi double layer"""
            res[offset + ptr['phi_dl']] = SV_dot[offset + ptr['phi_dl']] - (-i_Far + i_el_m - i_el_p)*cat.dyInv/cat.C_dl/A_C
            
            """Algebraic expression for charge neutrality in all phases"""
            res[offset + ptr['phi_ed']] = i_el_m - i_el_p + i_io_m - i_io_p
            
            """Calculate change in S8 nucleation sites"""
            res[offset + ptr['np_S8']] = SV_dot[offset + ptr['np_S8']]
            
            """Calculate change in Li2S nucleation sites"""
            res[offset + ptr['np_Li2S']] = SV_dot[offset + ptr['np_Li2S']]
        
        """==================Separator boundary conditions=================="""
        
        i_io_m = i_io_p
        N_io_m = N_io_p
        s1 = dict(s2)
        
        # Shift forward to NEXT node
        j = 0; offset = an.offsets[int(j)]
        s2 = set_state(SV, offset, an.ptr)
        
        dyInv_boundary = 1/(0.5*(sep.dy + an.dy))
        
        # Shift back to THIS node
        offset = sep.offsets[int(j)]
        
        D_el = sep.D_el 
        
        # Current node plus face boundary conditions
        N_io_p, i_io_p = dst(s1, s2, D_el, dyInv_boundary)
        
        res[offset + sep.ptr['rho_k_el']] = (SV_dot[offset + sep.ptr['rho_k_el']]
        - (N_io_m - N_io_p)*sep.dyInv)/sep.epsilon_el
        
        res[offset + sep.ptr['phi']] = i_io_m - i_io_p
        
        """====================Anode boundary conditions===================="""
          
        i_io_m = i_io_p
        N_io_m = N_io_p
        i_el_m = 0
        s1 = dict(s2)
        
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
        
        return res  
      
    "========================================================================="
    
    def set_state():
        return
    
    "========================================================================="
    
    def state_events(self, t, y, yd, sw):
        
        event1 = np.zeros([cat.npoints*int(y[cat.ptr_vec['np_S8']])])
        event2 = np.zeros([cat.npoints*int(y[cat.ptr_vec['np_S8']])])
        event1 = 1 - y[cat.ptr_vec['eps_S8']]
#        event2 = y[cat.ptr_vec['eps_S8']]
        
        event3 = np.zeros([cat.npoints*int(y[cat.ptr_vec['np_Li2S']])])
        event4 = np.zeros([cat.npoints*int(y[cat.ptr_vec['np_Li2S']])])
        event3 = 1 - y[cat.ptr_vec['eps_Li2S']]
        event4 = y[cat.ptr_vec['eps_Li2S']]
        
        event5 = np.zeros([cat.npoints])
#        event5 = 3.0 - y[cat.ptr_vec['phi_ed']]
        event6 = np.zeros([cat.npoints])
        event6 = y[cat.ptr_vec['phi_ed']] - 1.6
        
        
        events = np.concatenate((event1, event2, event3, event4, event5, event6))

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
    SV_eq, SV_dch, tags = main()
#    SV_eq_df, SV_ch_df, SV_req_df = main()
#    SV_eq, SV_ch, SV_req, SV_dch = main()

