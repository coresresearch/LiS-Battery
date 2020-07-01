# -*- coding: utf-8 -*-
"""
Created on Thu May  9 12:30:28 2019

@author: dkorff

This module contains the set up and running of the simulation for the Li-S
model.
"""

import numpy as np
import time
import cantera as ct
from matplotlib import pyplot as plt

# Import up the Assimulo DAE solver modules
from assimulo.solvers import IDA
from assimulo.problem import Implicit_Problem
from assimulo.exception import TerminateSimulation

# Initialize the model, parameters, and cantera objects:
from li_s_battery_init import anode, sep, cathode, inputs, sol_init, sim_time,\
    atol, rtol, sim_output, fig, axes

# Import post-processing routines:
from li_s_battery_post import label_columns, tag_strings, plot_sim, plot_mean_ps

def main():
    
    #TODO: Document what this line does, please.
    res_class = eval(inputs.test_type)
    
    # Save the current time, to measure the cpu time.
    t_count = time.time()
        
    "----------Equilibration----------"
    print('\nEquilibrating...')

    # Set external current to 0 for equilibration
    i_ext_eq = 0
    
    # Create problem object, simulation object, and run the simulation:
    SV_eq_df, SV_0, SV_dot_0 = run_model(res_class, sol_init.SV_0, 
        sol_init.SV_dot_0, i_ext_eq, sol_init.algvar, sim_time, atol, rtol, 
        sim_output)
    
    t_equilibrate = time.time() - t_count
    print("Equilibration time = ", t_equilibrate, '\n')
    
    print('Done equilibrating\n')
        
    for cycle_num in np.arange(0, inputs.n_cycles):
        "------------Discharging-------------"
        print('Discharging...')
        
        # Update problem object, simulation object, and run the simulation:
        SV_discharge_df, SV_0, SV_dot_0 = run_model(res_class, SV_0, SV_dot_0, 
            cathode.i_ext_amp, sol_init.algvar, sim_time, atol, rtol, 
            sim_output)
        
        print('Done Discharging\n')
        
        if inputs.flag_re_eq == 1:
            "--------Re-equilibration---------"
            print('Re-equilibrating...')
            # Set external current
            SV_re_eq_df, SV_0, SV_dot_0 = run_model(res_class, SV_0, SV_dot_0, 
                i_ext_eq, sol_init.algvar, sim_time, atol, rtol, sim_output)
        
            print('Done re-equilibrating\n')
            
        "-----------Charging-----------"
        print('Charging...')
        
        #TODO: Make this a conditional, that only sets the value if it exceeds 
        #   #the cutoff.
        SV_0[cathode.ptr_vec['eps_S8']] = cathode.eps_cutoff 

         # Set up and run the simulation:
        SV_charge_df, SV_0, SV_dot_0 = run_model(res_class, SV_0, SV_dot_0, 
                -cathode.i_ext_amp, sol_init.algvar, sim_time, atol, rtol, sim_output)
       
        print('Done Charging\n')    
    
    t_elapsed = time.time() - t_count
    print('t_cpu=', t_elapsed, '\n')    
    
    # Obtain tag strings for dataframe columns
    tags = tag_strings(SV_discharge_df)

    # Add results to the plot:
    plot_sim(tags, SV_discharge_df, 'Discharging', 0+2*cycle_num, fig, axes)
    plot_mean_ps(SV_discharge_df, tags, 'Discharging')

    if inputs.flag_re_eq:
        plot_sim(tags, SV_re_eq_df, 'Re-Equilibrating', 1, fig, axes)
        
    # Add results to the plots:
    plot_sim(tags, SV_charge_df, 'Charging', 
        1+inputs.flag_re_eq+2*cycle_num, fig, axes)
    plot_mean_ps(SV_charge_df, tags, 'Charging')

    return SV_eq_df, SV_discharge_df, SV_charge_df, tags 
    
"=============================================================================" 
"===========RESIDUAL CLASSES AND HELPER FUNCTIONS BEYOND THIS POINT==========="
"============================================================================="
from li_s_battery_init import sulfur_obj, Li2S_obj, carbon_obj, elyte_obj, \
    sulfur_elyte_surf_obj, Li2S_elyte_surf_obj, carbon_elyte_surf_obj, \
    Li2S_tpb_obj, lithium_obj, lithium_elyte_surf_obj, conductor_obj
from li_s_battery_functions import read_state_cathode, read_state_anode, \
    read_state_sep, set_geom, set_rxn, dst
from math import pi, exp, tanh

class cc_cycling(Implicit_Problem):    
    def res_fun(t, SV, SV_dot):
        
        res = np.zeros_like(SV)
        ptr = cathode.ptr; F = ct.faraday; R = ct.gas_constant; T = inputs.T
        
        """=============================CATHODE============================="""
        """CC BOUNDARY"""
        j = 0; offset = cathode.offsets[int(j)]
        i_ext = cathode.get_i_ext()
        s2, np_S8_2, np_Li2S_2, eps_S8_2, eps_Li2S_2, eps_elyte_2 = \
            read_state_cathode(SV, offset, cathode.ptr)

        # Set electronic current and ionic current boundary conditions. These 
        #   will be re-assigned as "in" values, once inside the for loop:
        i_el_out = i_ext
        i_io_out = 0
        N_k_elyte_out = np.zeros_like(elyte_obj.X)
        
        """=============================CATHODE============================="""
        """INTERIOR NODES"""
        for j in np.arange(1, cathode.npoints):
            
            # Set previous outlet fluxes to new inlet fluxes
            i_el_in = i_el_out
            i_io_in = i_io_out
            N_k_elyte_in = N_k_elyte_out

            # Set current state equal to previous "state 2" variables:
            s1 = dict(s2)
            np_S8 = np_S8_2; np_Li2S = np_Li2S_2; eps_S8 = eps_S8_2
            eps_Li2S = eps_Li2S_2; eps_elyte = eps_elyte_2
            
            # Update offset to NEXT node
            offset = cathode.offsets[int(j)]
            
            s2, np_S8_2, np_Li2S_2, eps_S8_2, eps_Li2S_2, eps_elyte_2 = \
                read_state_cathode(SV, offset, cathode.ptr)

            # Set variables to CURRENT NODE value
            offset = cathode.offsets[int(j-1)]
            
            # Calculate new particle radii based on new volume fractions
            A_S = 2*pi*np_S8*(3*eps_S8/2/np_S8/pi)**(2/3)
            A_L = 2*pi*np_Li2S*(3*eps_Li2S/2/np_Li2S/pi)**(2/3)
            
            r_S = 3*eps_S8/A_S
            r_L = 3*eps_Li2S/A_L
            
            tpb_len = 3*eps_Li2S/(r_L**2)
            
            A_C = cathode.A_C_0 - (pi*np_S8*r_S**2) - (pi*np_Li2S*r_L**2)

            # Set states for THIS node            
            carbon_obj.electric_potential = s1['phi_ed']
            elyte_obj.electric_potential = s1['phi_el'] 
            conductor_obj.electric_potential = s1['phi_ed']
            elyte_obj.X = s1['X_k']
            
            D_el = cathode.D_el*eps_elyte**(1.5)

            # Current node plus face boundary fluxes
            i_el_out = cathode.sigma_eff*(s1['phi_ed'] - s2['phi_ed'])*cathode.dyInv
            N_k_elyte_out, i_io_out = dst(s1, s2, D_el, cathode.dy, cathode.dy)
            
            sdot_C = carbon_elyte_surf_obj.get_net_production_rates(elyte_obj)
            R_C = sdot_C*A_C
            mult = tanh(eps_S8/cathode.eps_dropoff)  
            sdot_S8 = sulfur_elyte_surf_obj.get_creation_rates(sulfur_obj) - mult*sulfur_elyte_surf_obj.get_destruction_rates(sulfur_obj)
            sdot_S = sulfur_elyte_surf_obj.get_net_production_rates(elyte_obj)  
            R_S = sdot_S*A_S
    
            mult = tanh(eps_Li2S/cathode.eps_dropoff)  
            sdot_Li2S = Li2S_elyte_surf_obj.get_creation_rates(Li2S_obj) - mult*(Li2S_elyte_surf_obj.get_destruction_rates(Li2S_obj))
            sdot_L = Li2S_elyte_surf_obj.get_net_production_rates(elyte_obj)
            sdot_tpb = Li2S_tpb_obj.get_creation_rates(Li2S_obj) - mult*(Li2S_tpb_obj.get_destruction_rates(Li2S_obj))
            sdot_tpb_el = mult*Li2S_tpb_obj.get_creation_rates(elyte_obj) - Li2S_tpb_obj.get_destruction_rates(elyte_obj)
            R_L = sdot_L*A_L + sdot_tpb_el*tpb_len
    
            i_C = (carbon_elyte_surf_obj.get_net_production_rates(conductor_obj)*A_C + 
              (Li2S_tpb_obj.get_creation_rates(conductor_obj)*mult - 
               Li2S_tpb_obj.get_destruction_rates(conductor_obj))*tpb_len)
            i_Far = (i_C)*F/cathode.dyInv
            
            # Net rate of formation
            R_net = R_C + R_S + R_L 
            R_net[cathode.ptr['iFar']] += (-i_Far + i_el_in - i_el_out)/cathode.dy/F
            
            """Calculate change in Sulfur"""                
            res[offset + ptr['eps_S8']] = (SV_dot[offset + ptr['eps_S8']] 
                                        - sulfur_obj.volume_mole*sdot_S8*A_S)
       
            """Calculate change in Li2S"""
            res[offset + ptr['eps_Li2S']] = (SV_dot[offset + ptr['eps_Li2S']] 
                                          - Li2S_obj.volume_mole*(sdot_Li2S*A_L
                                          + sdot_tpb*tpb_len))
            
            """Calculate change in electrolyte"""
            res[offset + ptr['C_k_elyte']] = (SV_dot[offset + ptr['C_k_elyte']] - 
            (R_net + (N_k_elyte_in - N_k_elyte_out)*cathode.dyInv)/eps_elyte 
            + SV[offset + ptr['C_k_elyte']]*(- SV_dot[offset + ptr['eps_S8']] 
                                            - SV_dot[offset + ptr['eps_Li2S']])/eps_elyte)
            
            """Calculate change in delta-phi double layer"""
            res[offset + ptr['phi_dl']] = (SV_dot[offset + ptr['phi_dl']] - 
            (-i_Far + i_el_in - i_el_out)*cathode.dyInv/cathode.C_dl/A_C)
            
            """Algebraic expression for charge neutrality in all phases"""
            res[offset + ptr['phi_ed']] = i_el_in - i_el_out + i_io_in - i_io_out
            
            """Calculate change in S8 nucleation sites"""
            res[offset + ptr['np_S8']] = SV_dot[offset + ptr['np_S8']]
            
            """Calculate change in Li2S nucleation sites"""
            res[offset + ptr['np_Li2S']] = SV_dot[offset + ptr['np_Li2S']]
            
        """=============================CATHODE============================="""
        """SEPARATOR BOUNDARY"""
        
        i_el_in = i_el_out
        i_io_in = i_io_out
        N_k_elyte_in = N_k_elyte_out
        s1 = dict(s2)
        
        # Shift forward to NEXT node, first separator node (j=0)
        j = 0; offset = sep.offsets[int(j)]
        
        s2 = read_state_sep(SV, offset, sep.ptr)
        
        # Shift back to THIS node, set THIS node outlet conditions
        j = cathode.npoints-1; offset = cathode.offsets[int(j)]
        
        # Set variables to CURRENT NODE value
#        geom = set_geom(SV, offset, cathode.ptr)
        
        np_S8 = SV[offset + ptr['np_S8']]
        np_Li2S = SV[offset + ptr['np_Li2S']]
        eps_S8 = max(SV[offset + ptr['eps_S8']], cathode.eps_cutoff)
        eps_Li2S = max(SV[offset + ptr['eps_Li2S']], cathode.eps_cutoff)
        eps_elyte = 1 - cathode.eps_C_0 - eps_S8 - eps_Li2S  
        
        # Calculate new particle radii based on new volume fractions
        A_S = 2*pi*np_S8*(3*eps_S8/2/np_S8/pi)**(2/3)
        A_L = 2*pi*np_Li2S*(3*eps_Li2S/2/np_Li2S/pi)**(2/3)
        
        r_S = 3*eps_S8/A_S
        r_L = 3*eps_Li2S/A_L
        
        tpb_len = 3*eps_Li2S/(r_L**2)
        
        A_C = cathode.A_C_0 - (pi*np_S8*r_S**2) - (pi*np_Li2S*r_L**2)
        
        carbon_obj.electric_potential = s1['phi_ed']
        elyte_obj.electric_potential = s1['phi_el']
        conductor_obj.electric_potential = s1['phi_ed']
        elyte_obj.X = s1['X_k']
        
        # Set outlet boundary conditions for THIS node
        i_el_out = 0
        D_el = cathode.D_el*eps_elyte**(1.5)
        N_k_elyte_out, i_io_out = dst(s1, s2, D_el, cathode.dy, sep.dy)
        
        sdot_C = carbon_elyte_surf_obj.get_net_production_rates(elyte_obj)
        R_C = sdot_C*A_C
        mult = tanh(eps_S8/cathode.eps_dropoff)  
        sdot_S8 = sulfur_elyte_surf_obj.get_creation_rates(sulfur_obj) - mult*sulfur_elyte_surf_obj.get_destruction_rates(sulfur_obj)
        sdot_S = sulfur_elyte_surf_obj.get_net_production_rates(elyte_obj)  
        R_S = sdot_S*A_S

        mult = tanh(eps_Li2S/cathode.eps_dropoff)  
        sdot_Li2S = (Li2S_elyte_surf_obj.get_creation_rates(Li2S_obj) 
            - mult*(Li2S_elyte_surf_obj.get_destruction_rates(Li2S_obj)))
        sdot_L = Li2S_elyte_surf_obj.get_net_production_rates(elyte_obj)
        sdot_tpb = Li2S_tpb_obj.get_creation_rates(Li2S_obj) - mult*(Li2S_tpb_obj.get_destruction_rates(Li2S_obj))
        sdot_tpb_el = mult*Li2S_tpb_obj.get_creation_rates(elyte_obj) - Li2S_tpb_obj.get_destruction_rates(elyte_obj)
        R_L = sdot_L*A_L + sdot_tpb_el*tpb_len

        i_C = (carbon_elyte_surf_obj.get_net_production_rates(conductor_obj)*A_C + 
              (Li2S_tpb_obj.get_creation_rates(conductor_obj)*mult - 
               Li2S_tpb_obj.get_destruction_rates(conductor_obj))*tpb_len)
        i_Far = (i_C)*F/cathode.dyInv
        
        # Net rate of formation
        R_net = R_C + R_S + R_L 
        R_net[cathode.ptr['iFar']] += (-i_Far + i_el_in - i_el_out)/cathode.dy/F
                                 
        """Calculate change in Sulfur"""                
        res[offset + ptr['eps_S8']] = (SV_dot[offset + ptr['eps_S8']] 
                                    - sulfur_obj.volume_mole*sdot_S8*A_S)
        
        """Calculate change in Li2S"""
        res[offset + ptr['eps_Li2S']] = (SV_dot[offset + ptr['eps_Li2S']] 
                                      - Li2S_obj.volume_mole*(sdot_Li2S*A_L
                                      + sdot_tpb*tpb_len))
                
        """Calculate change in electrolyte"""
        res[offset + ptr['C_k_elyte']] = (SV_dot[offset + ptr['C_k_elyte']] - 
        (R_net + (N_k_elyte_in - N_k_elyte_out)*cathode.dyInv)/eps_elyte
        + SV[offset + ptr['C_k_elyte']]*(- SV_dot[offset + ptr['eps_S8']] 
                                        - SV_dot[offset + ptr['eps_Li2S']])/eps_elyte)
        
        """Calculate change in delta-phi double layer"""
        res[offset + ptr['phi_dl']] = (SV_dot[offset + ptr['phi_dl']] - 
        (-i_Far + i_el_in - i_el_out)*cathode.dyInv/cathode.C_dl/A_C)
        
        """Algebraic expression for charge neutrality in all phases"""
        res[offset + ptr['phi_ed']] = i_el_in - i_el_out + i_io_in - i_io_out
        
        """Calculate change in S8 nucleation sites"""
        res[offset + ptr['np_S8']] = SV_dot[offset + ptr['np_S8']]
        
        """Calculate change in Li2S nucleation sites"""
        res[offset + ptr['np_Li2S']] = SV_dot[offset + ptr['np_Li2S']]
        
        """============================SEPARATOR============================"""
        """INTERIOR NODES"""
        
        """============================SEPARATOR============================"""
        """CATHODE BOUNDARY"""
        
        i_io_in = i_io_out
        N_k_elyte_in = N_k_elyte_out
        s1 = dict(s2)
        
        # Shift forward to NEXT node
        j = 0; offset = anode.offsets[int(j)]
        s2 = read_state_anode(SV, offset, j, anode.ptr)
                
        # Shift back to THIS node
        offset = sep.offsets[int(j)]
        
        D_el = sep.D_el 
        
        # Current node plus face boundary conditions
        N_k_elyte_out, i_io_out = dst(s1, s2, D_el, sep.dy, anode.dy)
        
        res[offset + sep.ptr['C_k_elyte']] = (SV_dot[offset + sep.ptr['C_k_elyte']]
        - (N_k_elyte_in - N_k_elyte_out)*sep.dyInv/sep.epsilon_el)
                
        res[offset + sep.ptr['phi']] = i_io_in - i_io_out
        
        """==============================ANODE=============================="""
        """INTERIOR NODES"""
          
        i_io_in = i_io_out
        N_k_elyte_in = N_k_elyte_out
        i_el_in = 0
        s1 = dict(s2)
        
        j = 0
        offset = anode.offsets[int(j)]
        
        i_el_out = i_ext
        i_io_out = 0
        N_k_elyte_out = 0
        
        elyte_obj.X = s1['X_k']
        elyte_obj.electric_potential = s1['phi_el']
        lithium_obj.electric_potential = s1['phi_ed']
        conductor_obj.electric_potential = s1['phi_ed']
        
        sdot_Li = lithium_elyte_surf_obj.get_net_production_rates(elyte_obj)
        sdot_Far = lithium_elyte_surf_obj.get_net_production_rates(conductor_obj)
        
        R_net = sdot_Li*anode.A_Li
        i_Far = sdot_Far*anode.A_Li*F*anode.dy
                   
        res[anode.ptr['C_k_elyte'][j]] = (SV_dot[anode.ptr['C_k_elyte'][j]]
            - (R_net + (N_k_elyte_in - N_k_elyte_out)*anode.dyInv)/anode.eps_el)

        res[anode.ptr['phi_dl'][j]] = (
            SV_dot[anode.ptr['phi_dl'][j]]
            - (-i_Far + i_el_in - i_el_out)*anode.dyInv/anode.C_dl/anode.A_Li) 
        
        res[anode.ptr['phi_ed'][j]] = SV[anode.ptr['phi_ed'][j]]
        
        """==============================ANODE=============================="""
        """CC BOUNDARY"""

        
        return res  
    
    "========================================================================="
    
    def state_events(self, t, y, yd, sw):
        
        event1 = np.zeros([cathode.npoints])
        event2 = np.zeros([cathode.npoints])
        event1 = 1 - y[cathode.ptr_vec['eps_S8']]
#        event2 = y[cathode.ptr_vec['eps_S8']]
        
        event3 = np.zeros([cathode.npoints])
        event4 = np.zeros([cathode.npoints])
        event3 = 1 - y[cathode.ptr_vec['eps_Li2S']]
#        event4 = y[cathode.ptr_vec['eps_Li2S']]
        
        event5 = np.zeros([cathode.npoints])
        event5 = 2.8 - y[cathode.ptr_vec['phi_ed']]
        event6 = np.zeros([cathode.npoints])
        event6 = y[cathode.ptr_vec['phi_ed']] - 1.5
        
        event7 = np.zeros([cathode.npoints*elyte_obj.n_species])
        event7 = y[cathode.ptr_vec['C_k_elyte']]
        
        
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
            raise TerminateSimulation
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
    
def create_problem(res_class, SV_0, SV_dot_0, algvar, t_0):
    bat = res_class(res_class.res_fun, SV_0, SV_dot_0, t_0)
    bat.external_event_detection = True
    bat.algvar = algvar

    return bat

def run_model(res_class, SV_0, SV_dot_0, i_ext, algvar, time, 
    atol, rtol, output):

    # Set external current to 0 for equilibration
    cathode.set_i_ext(i_ext)
    
    # Create problem object and simulation object, then run the simulation:
    bat = res_class(res_class.res_fun, SV_0, SV_dot_0, time[0])
    bat.external_event_detection = True
    bat.algvar = algvar
    
    sim = IDA(bat)
    sim.atol = atol
    sim.rtol = rtol
    sim.verbosity = output
    sim.make_consistent('IDA_YA_YDP_INIT')

    t, SV, SV_dot = sim.simulate(time[1])
    
    SV_0 = SV[-1, :]
    SV_dot_0 = SV_dot[-1, :]

    # Put solution into pandas dataframe with labeled columns
    SV_df = label_columns(t, SV, anode.npoints, sep.npoints, cathode.npoints)

    return SV_df, SV_0, SV_dot_0

def create_sim(bat, atol, rtol, output):
    sim = IDA(bat)
    sim.atol = atol
    sim.rtol = rtol
    sim.verbosity = output
    sim.make_consistent('IDA_YA_YDP_INIT')

    return sim

if __name__ == "__main__":
    SV_eq, SV_discharge, SV_ch, tags = main()

"============================================================================="
