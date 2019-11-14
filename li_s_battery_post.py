# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 14:35:17 2019

@author: dkorff
"""
from li_s_battery_inputs import inputs
from li_s_battery_init import anode
from li_s_battery_init import cathode
from li_s_battery_init import sep
from li_s_battery_init import elyte_obj, sulfur_obj, Li2S_obj, carbon_obj
from matplotlib import pyplot as plt

import numpy as np
import pandas as pd

def plot_sim(tags, SV_df, stage, yax, fig, axes):
    
    if stage == 'Discharging':
        showlegend = 1
    else:
        showlegend = 0
    
#    SV_df = SV_df_orig.copy()
#    SV_df['phi_dl'] = SV_df['phi_dl'] + SV_df['phi_el']
    
    vol_fracs = tags['eps_S8'] + tags['eps_Li2S']
#    phi = tags['phi_dl'] + tags['phi_ed']
    phi = tags['phi_ed']
    fontsize = 18
   
    SV_df.loc[:, 'Time'] *= -cathode.i_ext_amp*inputs.A_cat/3600/cathode.m_S_0
    t = SV_df['Time']
    # Plot potential for the electrolyte and the double layer
    SV_plot = SV_df.plot(x='Time', y=phi, ax=axes[0],  xlim=[0,t.iloc[-1]])
    SV_plot.set_title(stage, fontsize = fontsize)
    SV_plot.set_ylabel(r'$V_{cell}$ [V]', fontsize = fontsize)
    SV_plot.set_xlabel('Capacity $[Ah/kg_{sulfur}]$', fontsize = fontsize).set_visible(False)
    SV_plot.set_xlim((0, 1750))
#    SV_plot.set_ylim((2.25, 2.5))
    SV_plot.legend(loc=2, bbox_to_anchor=(1.0, 1), ncol=1, borderaxespad=0,
                   frameon=False, fontsize = 15).set_visible(False)
    SV_plot.tick_params(axis='both', labelsize=16)
#    SV_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
    # Plot Li2S and S8 volume fractions
    SV_plot = SV_df.plot(x='Time', y=vol_fracs, ax=axes[1], xlim=[0,t.iloc[-1]])
#    SV_plot.set_title(stage, fontsize = fontsize)
    SV_plot.set_ylabel(r'$\varepsilon_i$ [-]', fontsize = fontsize)
    SV_plot.set_xlabel('Time [s]', fontsize = fontsize).set_visible(False)
    SV_plot.set_xlim((0, 1750))
    SV_plot.legend(loc=2, bbox_to_anchor=(1.0, 1), ncol=1, borderaxespad=0,
                   frameon=False, fontsize = 15).set_visible(showlegend)
    SV_plot.tick_params(axis='both', labelsize=16)
#    SV_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
    # Plot species densities in electrolyte
    SV_plot = SV_df.plot(x='Time', y=tags['rho_el'][4:], logy=True, ax=axes[2], xlim=[0,t.iloc[-1]]) #
#    SV_plot.set_title(stage, fontsize = fontsize)
    SV_plot.set_ylabel(r'$\rho_k$ [kmol/m$^3]$', fontsize = fontsize)
    SV_plot.set_xlabel('Capacity $[Ah/kg_{sulfur}]$', fontsize = fontsize).set_visible(True)
    SV_plot.set_ylim((1e-12, 1e2))
    SV_plot.set_xlim((0, 1750))
    SV_plot.legend(loc=2, bbox_to_anchor=(1.0, 1), ncol=1, borderaxespad=0,
                   frameon=False, fontsize = 15).set_visible(showlegend)
    SV_plot.tick_params(axis='both', labelsize=16)
#    SV_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
#    # Plot species densities in separator electrolyte
#    SV_plot = SV_df.plot(x='Time', y=tags['rho_el_sep'][4:], logy=True, ax=axes[3], xlim=[0,t.iloc[-1]]) #
##    SV_plot.set_title(stage, fontsize = fontsize)
#    SV_plot.set_ylabel(r'$\rho_k$ [kmol/m$^3]$', fontsize = fontsize)
#    SV_plot.set_xlabel('Capacity $[Ah/kg_{sulfur}]$', fontsize = fontsize).set_visible(True)
#    SV_plot.legend(loc=2, bbox_to_anchor=(1.0, 1), ncol=1, borderaxespad=0,
#                   frameon=False, fontsize = 15).set_visible(showlegend)
#    SV_plot.tick_params(axis='both', labelsize=16)
##    SV_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#    
#    # Plot species densities in anode electrolyte
#    SV_plot = SV_df.plot(x='Time', y=tags['rho_el_an'][4:], logy=True, ax=axes[4], xlim=[0,t.iloc[-1]]) #
##    SV_plot.set_title(stage, fontsize = fontsize)
#    SV_plot.set_ylabel(r'$\rho_k$ [kmol/m$^3]$', fontsize = fontsize)
#    SV_plot.set_xlabel('Capacity $[Ah/kg_{sulfur}]$', fontsize = fontsize).set_visible(True)
#    SV_plot.legend(loc=2, bbox_to_anchor=(1.0, 1), ncol=1, borderaxespad=0,
#                   frameon=False, fontsize = 15).set_visible(showlegend)
#    SV_plot.tick_params(axis='both', labelsize=16)
##    SV_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
    
    return

"============================================================================="

def label_columns(t, SV, an_np, sep_np, cat_np):
    
    # Convert t and SV arrays into pandas data frames
    t_df = pd.DataFrame(t)
    SV_df = pd.DataFrame(SV)
    
    # Set the column label for the t data frame to the number of columns in SV
    newcols_t = {0: SV_df.shape[1]}
    t_df.rename(columns = newcols_t, inplace = True)
    
    # Concatenate t_df onto end of SV_df by columns (axis = 1)
    SV_df = pd.concat((SV_df, t_df), axis = 1)
    
    """Label anode points"""
    newcols = {}
    for j in np.arange(0, an_np):
        offset = anode.offsets[j]  # Set node offset value for loop
        
#        # Loop over number of shells in anode
#        for k in np.arange(0, anode.nshells):
#            newcols_an = {k + offset: 'X_an'+str(j+1)+str(k+1)}
#            newcols.update(newcols_an)
            
        # Loop over number of species in electrolyte
        for k in np.arange(0, elyte_obj.n_species):
            species = elyte_obj.species_names[k]
            newcols_el = {k + offset: 'X_'+species+'_an'+str(j+1)}
            newcols.update(newcols_el)
            
        # Add tags for electrod and double layer potentials
        newcols_phi = {0+elyte_obj.n_species+offset: 'Phi_an_dl'+str(j+1),
                       1+elyte_obj.n_species+offset: 'Phi_an'+str(j+1)}
        newcols.update(newcols_phi)
        
        SV_df.rename(columns=newcols, inplace = True)
        
        """Label separator points"""
        newcols = {}
    for j in np.arange(0, sep_np):
        offset = sep.offsets[j] # Set node offset value for loop
        
        # Loop over number of species in electrolyte
        for k in np.arange(0, elyte_obj.n_species):
            species = elyte_obj.species_names[k]
            newcols_el = {k + offset: 'X_'+species+'_sep'+str(j+1)}
            newcols.update(newcols_el)
            
        # Add tag for electrolyte potential
        newcols_phi = {0+elyte_obj.n_species+offset: 'Phi_sep'+str(j+1)}
        newcols.update(newcols_phi)
        
        SV_df.rename(columns=newcols, inplace = True)
        
    """Label cathode points"""
    newcols = {}
    for j in np.arange(0, cat_np):
        offset = cathode.offsets[j]  # Set node offset value for loop
        
        # Add tags for particle radius of Li2S and S8
        newcols_r = {0+offset: 'eps_S8'+str(j+1),
                     1+offset: 'eps_Li2S'  +str(j+1)}
        newcols.update(newcols_r)
        
        # Loop over number of species in electrolyte
        for k in np.arange(0, elyte_obj.n_species):
            spec = elyte_obj.species_names[k]
            newcols_el = {2 + k + offset: 'rho_'+spec+'_cat'+str(j+1)}
            newcols.update(newcols_el)
            
        # Add tags for double layer and electrolyte potentials
        newcols_phi = {2 + elyte_obj.n_species + offset: 'Phi_dl'+str(j+1),
                       3 + elyte_obj.n_species + offset: 'Phi_ed'+str(j+1)}
        newcols.update(newcols_phi)
        
        SV_df.rename(columns = newcols, inplace = True)
        
        # Add tag for number of nucleation sites
        newcols_nucl = {4 + elyte_obj.n_species + offset: 'np_S8'+str(j+1),
                        5 + elyte_obj.n_species + offset: 'np_Li2S'+str(j+1)}
        newcols.update(newcols_nucl)
        
        SV_df.rename(columns = newcols, inplace = True)
        
    newcols_time = {SV_df.shape[1]-1: 'Time'}
    SV_df.rename(columns = newcols_time, inplace = True)
    
    return SV_df

"============================================================================="

def tag_strings(SV):
    
    SV_labels = SV.columns.values.tolist()
    
    r_Li2S = np.array([])
    r_S8 = np.array([])
    rho_el = []
    phi_dl = np.array([])
    phi_ed = np.array([])
    np_S8 = np.array([])
    np_Li2S = np.array([])
    
    rho_el_sep = []
    phi_sep = np.array([])
    
    rho_el_an = []
    phi_dl_an = np.array([])
    phi_an = np.array([])
    
    ptr = cathode.ptr
    for j in np.arange(0, cathode.npoints):
        offset = int(cathode.offsets[j])
        
        r_Li2S = np.append(r_Li2S, SV_labels[ptr['eps_Li2S'] + offset])
        r_S8 = np.append(r_S8, SV_labels[ptr['eps_S8'] + offset])
        
        rho_el[0 + offset:elyte_obj.n_species + offset] = \
            SV_labels[ptr['rho_k_el'][0]+offset:ptr['rho_k_el'][-1]+offset+1]
            
        phi_dl = np.append(phi_dl, SV_labels[ptr['phi_dl'] + offset])
        phi_ed = np.append(phi_ed, SV_labels[ptr['phi_ed'] + offset])
        np_S8 = np.append(np_S8, SV_labels[ptr['np_S8'] + offset])
        np_Li2S = np.append(np_Li2S, SV_labels[ptr['np_Li2S'] + offset])
        
    ptr = sep.ptr
    for j in np.arange(0, sep.npoints):
        offset = int(sep.offsets[j])
        
        rho_el_sep[0 + offset:elyte_obj.n_species + offset] = \
            SV_labels[ptr['rho_k_el'][0]+offset:ptr['rho_k_el'][-1]+offset+1]
            
        phi_sep = np.append(phi_sep, SV_labels[ptr['phi'] + offset])
        
    ptr = anode.ptr
    for j in np.arange(0, anode.npoints):
        offset = int(anode.offsets[j])
        
        rho_el_an[0 + offset:elyte_obj.n_species + offset] = \
            SV_labels[ptr['rho_k_el'][0]+offset:ptr['rho_k_el'][-1]+offset+1]
            
        phi_dl_an = np.append(phi_dl_an, SV_labels[ptr['phi_dl'] + offset])
        phi_an = np.append(phi_an, SV_labels[ptr['phi_ed'] + offset])
        
    phi_sep = phi_sep.tolist()
    phi_dl_an = phi_dl_an.tolist()
    phi_an = phi_an.tolist()
        
    r_Li2S = r_Li2S.tolist()
    r_S8 = r_S8.tolist()
    phi_dl = phi_dl.tolist()
    phi_ed = phi_ed.tolist()
    np_S8 = np_S8.tolist()
    np_Li2S = np_Li2S.tolist()
    
    tags = {}
    tags['eps_Li2S'] = r_Li2S; tags['eps_S8'] = r_S8; tags['rho_el'] = rho_el
    tags['phi_dl'] = phi_dl; tags['phi_ed'] = phi_ed; tags['np_S8'] = np_S8
    tags['np_Li2S'] = np_Li2S; tags['rho_el_sep'] = rho_el_sep; tags['phi_sep'] = phi_sep
    tags['rho_el_an'] = rho_el_an; tags['phi_dl_an'] = phi_dl_an; tags['phi_an'] = phi_an
    
    return tags
    
    
    
    
    
    
    
    
    
    
    
    