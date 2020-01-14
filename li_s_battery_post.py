# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 14:35:17 2019

@author: dkorff
"""
from li_s_battery_inputs import inputs
from li_s_battery_init import anode
from li_s_battery_init import cathode
from li_s_battery_init import sep
from li_s_battery_init import elyte_obj, sulfur_obj, Li2S_obj, carbon_obj, conductor_obj
from li_s_battery_init import carbon_el_s, Li2S_el_s, sulfur_el_s
from li_s_battery_functions import dst, set_state, set_state_sep
from matplotlib import pyplot as plt
from math import pi
import numpy as np
import pandas as pd
import cantera as ct

def conservation_tests(SV, tags):
    F = ct.faraday
    flag_cat = 1
    flag_sep = 1
    flag_an = 1
    
    n_S_0 = flag_cat*cathode.n_S_0 + flag_sep*sep.n_S_0 + flag_an*anode.n_S_0
    i_sep = np.zeros([len(SV.index)])
    N_Li_sep = np.zeros([len(SV.index)])
    i_cc = np.zeros([len(SV.index)])
    i_dl = np.zeros([len(SV.index)])
    n_S_tot = np.zeros([len(SV.index)])
    n_S_elyte = np.zeros([len(SV.index)])
    n_S_solid_vec = np.zeros([len(SV.index)])
    n_S_Li2S_vec = np.zeros([len(SV.index)])
    n_S_cat = np.zeros([len(SV.index)])
    n_S_sep = np.zeros([len(SV.index)])
    n_S_an = np.zeros([len(SV.index)])
    charge_el_cat = np.zeros([len(SV.index)])
    charge_el_sep = np.zeros([len(SV.index)])
    charge_el_an = np.zeros([len(SV.index)])
    i_ext = cathode.i_ext_amp
    n_Li_cat = np.zeros([len(SV.index)])
    n_Li_Li2S = np.zeros([len(SV.index)])
    n_Li_tot = np.zeros([len(SV.index)])
    N_Li_dl = np.zeros([len(SV.index)])
    N_Li_dl_integral = np.zeros([len(SV.index)])
    
    eps_S_vec = np.zeros([len(SV.index)])
    eps_Li2S_vec = np.zeros([len(SV.index)])
    eps_el_vec = np.zeros([len(SV.index)])
    eps_C_vec = np.zeros([len(SV.index)])
    
    for i, state in SV.iterrows():
        # We will check several items to ensure conservation at each time step.
        # All quantities will be per cell area
        #   1. sulfur atoms in all phases
        #   2. lithium atoms in all phases
        #   3. charge neutrality in electrolyte and carbon
        
        """1. Conservation of sulfur"""
        np_S = state.iloc[cathode.ptr['np_S8']]
        np_L = state.iloc[cathode.ptr['np_Li2S']]
        
        # Volume fractions at current state
        eps_S8 = max(state.iloc[cathode.ptr['eps_S8']], 1e-25)
        eps_Li2S = max(state.iloc[cathode.ptr['eps_Li2S']], 1e-25)
        eps_el = 1 - eps_S8 - eps_Li2S - cathode.eps_C_0
        
        eps_S_vec[i] = eps_S8
        eps_Li2S_vec[i] = eps_Li2S
        eps_el_vec[i] = eps_el
        eps_C_vec[i] = 1 - eps_S8 - eps_Li2S - eps_el

        # Concentration vector for all species in elyte at current state
        rho_el_cat = state.iloc[cathode.ptr['rho_k_el']]
        rho_el_sep = state.iloc[sep.offsets[int(0)]+sep.ptr['rho_k_el']]
        rho_el_an = state.iloc[anode.offsets[int(0)]+anode.ptr['rho_k_el']]
        
        # Concentration of just sulfur containing species in electrolyte
        rho_S_el_cat = state.iloc[cathode.ptr['rho_k_el'][cathode.i_S8:]]
        rho_S_el_sep = state.iloc[sep.offsets[int(0)]+sep.ptr['rho_k_el'][cathode.i_S8:]]
        rho_S_el_an = state.iloc[anode.offsets[int(0)]+anode.ptr['rho_k_el'][cathode.i_S8:]]

        # Number of moles of sulfur atoms in elyte of each component
        n_S_cat[i] = eps_el*cathode.H*np.dot(cathode.n_S_atoms, rho_S_el_cat)
        n_S_sep[i] = sep.epsilon_el*sep.H*np.dot(cathode.n_S_atoms, rho_S_el_sep)
        n_S_an[i] = anode.eps_el*anode.H*np.dot(cathode.n_S_atoms, rho_S_el_an)
        
        # Number of moles of sulfur atoms in solid phases
        n_S_solid = 8*sulfur_obj.density_mole*eps_S8*cathode.H
        n_S_Li2S = Li2S_obj.density_mole*eps_Li2S*cathode.H
        
        n_S_elyte[i] = n_S_cat[i] + n_S_sep[i] + n_S_an[i]
        n_S_solid_vec[i] = n_S_solid
        n_S_Li2S_vec[i] = n_S_Li2S
        n_S_tot[i] = (n_S_cat[i] + n_S_solid + n_S_Li2S) + n_S_sep[i] + n_S_an[i]
        
        """2. Conservation of lithium"""
        offset1 = cathode.offsets[-1]
        s1 = set_state(state.values, offset1, cathode.ptr)
        offset2 = sep.offsets[0]
        s2 = set_state_sep(state.values, offset2, sep.ptr)
        dyInv = 1/(0.5*(cathode.dy + sep.dy))
        D_el = cathode.D_el*eps_el**(1.5)
        N_io_sep, i_io_sep = dst(s1, s2, D_el, dyInv)
#        print(N_io_sep, '\n')
        if i == 0:
            dt = 0
        else:
            dt = SV.iloc[i, -1] - SV.iloc[i-1, -1]
            
        """Cantera objects (to get faradaic current)"""
        carbon_obj.electric_potential = s1['phi_ed']
        elyte_obj.electric_potential = s1['phi_el'] 
        conductor_obj.electric_potential = s1['phi_ed']
        
        elyte_obj.X = s1['X_k']
        
        # Calculate new particle radii based on new volume fractions
        A_S = 3*eps_S8/(3*eps_S8*cathode.V_0/2/pi/np_S)**(1/3)
        A_L = 3*eps_Li2S/(3*eps_Li2S*cathode.V_0/2/pi/np_L)**(1/3)
        
        r_S = 3*eps_S8/A_S
        r_L = 3*eps_Li2S/A_L
        
        A_C = inputs.A_C_0 - (pi*np_S*r_S**2)/cathode.V_0 - (pi*np_L*r_L**2)/cathode.V_0
        
        i_Far = carbon_el_s.get_net_production_rates(conductor_obj)*F*A_C/cathode.dyInv
            
        # Net rate of formation
        R_Li_dl = (-i_Far + i_ext - 0)/cathode.H/F
#        print(R_Li_dl)
        N_Li_dl[i] = dt*R_Li_dl*cathode.H

        N_Li_sep[i] = N_io_sep[2]*dt
        i_sep[i] = i_io_sep
        i_cc[i] = i_ext
        i_dl[i] = (-i_Far + i_ext - 0)
        
        rho_Li_el_cat = rho_el_cat[2]
        n_Li_cat[i] = eps_el*cathode.H*(rho_Li_el_cat - inputs.C_k_el_0[2])
        
        n_Li_Li2S[i] = 2*Li2S_obj.density_mole*(eps_Li2S - cathode.eps_L_0)*cathode.H
        
        n_Li_tot[i] = n_Li_cat[i] + n_Li_Li2S[i]
        
        """3. Charge neutrality"""
        charge_el_cat[i] = eps_el*cathode.H*np.dot(inputs.z_k_el, rho_el_cat)
        charge_el_sep[i] = sep.epsilon_el*sep.H*np.dot(inputs.z_k_el, rho_el_sep)
        charge_el_an[i] = anode.eps_el*anode.H*np.dot(inputs.z_k_el, rho_el_an)
#            C_S_anions_0 = inputs.C_k_el_0[5:]
#            C_S_anions = SV[offset + ptr['rho_k_el'][5:]]
#            
#            Q = -i_ext*t/F
#            Q_S = sum((C_S_anions - C_S_anions_0)*(-2)*eps_el*cat.H)
#            Q_dl = cat.C_dl*A_C*cat.H*SV[offset + ptr['phi_dl']]/F
#            print(Q, Q_S, Q_dl, Q + Q_S - Q_dl, t, '\n')
        
    pct_error_S = 100*(n_S_tot - n_S_0)/n_S_0
    N_Li_integral = np.cumsum(N_Li_sep) #*dt_vec
    N_Li_dl_integral = np.cumsum(N_Li_dl)
    pct_error_Li = (N_Li_integral + n_Li_tot)
    
    test = (n_S_0 - n_S_tot[-1])
    """---------------------------------------------------------------------"""
    """Plotting"""
    
    "-----Plot conservation of sulfur-----"
    fig = plt.figure(1)
    ax = fig.add_axes([0.2,0.2,0.6,0.75])
    fig.set_size_inches((10.,5.0))
    
    "Formatting for the figure:"
    fs = 20     #font size for plots
    lw = 2.0    #line width for plots
#    font = plt.matplotlib.font_manager.FontProperties(family='Times New Roman',size=fs-1)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')    
    
    p1, = plt.plot(SV.loc[:, 'Time'], n_S_tot, 'k-', linewidth=lw)
#    p2, = plt.plot(SV.loc[:, 'Time'], n_S_elyte, 'g-', linewidth=lw)
    p3, = plt.plot(SV.loc[:, 'Time'], n_S_solid_vec, linewidth=lw)
    p4, = plt.plot(SV.loc[:, 'Time'], n_S_Li2S_vec, linewidth=lw)
    p5, = plt.plot(SV.loc[:, 'Time'], n_S_0*np.ones((len(n_S_tot))), 'k--', linewidth=lw)
    p6, = plt.plot(SV.loc[:, 'Time'], n_S_cat, linewidth=lw)
    p7, = plt.plot(SV.loc[:, 'Time'], n_S_sep, linewidth=lw)
    p8, = plt.plot(SV.loc[:, 'Time'], n_S_an, linewidth=lw)
    plt.legend(['Total', 'S8', 'Li2S', 'Initial', 'Cathode elyte', 'Sep elyte', 'Anode elyte'])
#    p1, = plt.plot(SV_df.loc[:, 'Time'], SV_df.loc[:, tags['phi_ed']], 'k-', linewidth=lw)
#    plt.xlim((0, SV.loc[-1, 'Time']))
    plt.xticks([0, 30000, 60000, 90000, 120000, 150000, 180000])
#    plt.ylim((1.6, 2.6))
#    plt.yticks([2, 3, 4, 5, 6, 7, 8])
    plt.ylabel(r'$C_{\mathrm{S}}^{''} \hspace{0.5} [\mathrm{kmol}_{\mathrm{S}} \hspace{0.5} \mathrm{m}^{-2}]$', 
               fontstyle='normal', fontname='Times new Roman', fontsize=fs+2, labelpad=5.0)
#    plt.xlabel(r'Capacity $[\mathrm{Ah} \hspace{0.5} \mathrm{kg}^{-1}_{\mathrm{sulfur}}]$', \
#                fontstyle='normal', fontname='Times new Roman', fontsize=fs+2, labelpad=5.0)
    
    "-----Plot percent error in sulfur-----"
    fig = plt.figure(2)
    ax = fig.add_axes([0.2,0.2,0.6,0.75])
    fig.set_size_inches((10.,5.0))
    
    "Formatting for the figure:"
    fs = 20     #font size for plots
    lw = 2.0    #line width for plots
#    font = plt.matplotlib.font_manager.FontProperties(family='Times New Roman',size=fs-1)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')    
    
    p1, = plt.plot(SV.loc[:, 'Time'], pct_error_S, 'k-', linewidth=lw)
    p2, = plt.plot(SV.loc[:, 'Time'], np.zeros_like(pct_error_S), 'k--', linewidth=lw)
#    p1, = plt.plot(SV_df.loc[:, 'Time'], SV_df.loc[:, tags['phi_ed']], 'k-', linewidth=lw)
#    plt.xlim((0, SV.loc[-1, 'Time']))
#    plt.xticks([0, 250, 500, 750, 1000, 1250, 1500, 1750])
    plt.ylim((-10, 100))
    plt.xticks([0, 30000, 60000, 90000, 120000, 150000, 180000])
#    plt.yticks([2, 3, 4, 5, 6, 7, 8])
#    plt.ylabel('Mean PS order', fontstyle='normal', fontname='Times new Roman', \
#                fontsize=fs+2, labelpad=5.0)
#    plt.xlabel(r'Capacity $[\mathrm{Ah} \hspace{0.5} \mathrm{kg}^{-1}_{\mathrm{sulfur}}]$', \
#                fontstyle='normal', fontname='Times new Roman', fontsize=fs+2, labelpad=5.0)
    
    "-----Plot charge neutrality-----"
    fig = plt.figure(3)
    ax = fig.add_subplot(311)
    ax2 = fig.add_subplot(312, sharex=ax, sharey=ax)
    ax3 = fig.add_subplot(313, sharex=ax, sharey=ax)
    fig.set_size_inches((10.,8.0))
    
    "Formatting for the figure:"
    fs = 20     #font size for plots
    lw = 2.0    #line width for plots
#    font = plt.matplotlib.font_manager.FontProperties(family='Times New Roman',size=fs-1)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')    
    
    plt.sca(ax)
    p1, = plt.plot(SV.loc[:, 'Time'], charge_el_cat, 'k-', linewidth=lw)
    p2, = plt.plot(SV.loc[:, 'Time'], np.zeros_like(charge_el_cat), 'k--', linewidth=lw)
    plt.sca(ax2)
    p3, = plt.plot(SV.loc[:, 'Time'], charge_el_sep, 'k-', linewidth=lw)
    p4, = plt.plot(SV.loc[:, 'Time'], np.zeros_like(charge_el_sep), 'k--', linewidth=lw)
    plt.sca(ax3)
    p5, = plt.plot(SV.loc[:, 'Time'], charge_el_an, 'k-', linewidth=lw)
    p6, = plt.plot(SV.loc[:, 'Time'], np.zeros_like(charge_el_an), 'k--', linewidth=lw)
#    p1, = plt.plot(SV_df.loc[:, 'Time'], SV_df.loc[:, tags['phi_ed']], 'k-', linewidth=lw)
#    plt.xlim((0, SV.loc[-1, 'Time']))
#    plt.xticks([0, 250, 500, 750, 1000, 1250, 1500, 1750])
    plt.xticks([0, 30000, 60000, 90000, 120000, 150000, 180000])
#    plt.ylim((0, 100))
#    plt.yticks([2, 3, 4, 5, 6, 7, 8])
#    plt.ylabel('Mean PS order', fontstyle='normal', fontname='Times new Roman', \
#                fontsize=fs+2, labelpad=5.0)
#    plt.xlabel(r'Capacity $[\mathrm{Ah} \hspace{0.5} \mathrm{kg}^{-1}_{\mathrm{sulfur}}]$', \
#                fontstyle='normal', fontname='Times new Roman', fontsize=fs+2, labelpad=5.0)
    
    "-----Plot lithium balance in cathode-----"
    fig=plt.figure(4)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.75])
    fig.set_size_inches((10., 5.0))
    
    fs = 20
    lw = 2.0
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')    
        
    p1, = plt.plot(SV.loc[:, 'Time'], -i_sep, 'k-', linewidth=lw)
    p2, = plt.plot(SV.loc[:, 'Time'], i_cc, 'g-', linewidth=lw)
    p3, = plt.plot(SV.loc[:, 'Time'], i_cc - i_sep, linewidth=lw)
    plt.legend(['Separator', 'External', 'Net'])
    
    "-----Plot lithium balance in cathode-----"
    fig=plt.figure(5)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.75])
    fig.set_size_inches((10., 5.0))
    
    fs = 20
    lw = 2.0
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')    
        
    N_Li_flux_in = -N_Li_integral + N_Li_dl_integral
    p1, = plt.plot(SV.loc[:, 'Time'], N_Li_flux_in, 'k-', linewidth=lw)
    p2, = plt.plot(SV.loc[:, 'Time'], n_Li_cat, 'g-', linewidth=lw)
    p3, = plt.plot(SV.loc[:, 'Time'], n_Li_Li2S, linewidth=lw)
    p4, = plt.plot(SV.loc[:, 'Time'], n_Li_tot, linewidth=lw)
#    p5, = plt.plot(SV.loc[:, 'Time'], N_Li_dl_integral, linewidth=lw)
    plt.legend(['Lithium from separator and double layer', 'Change in elyte', 'Change in solid', 'Total change'])
    plt.xticks([0, 30000, 60000, 90000, 120000, 150000, 180000])
    
    "-----Plot lithium balance in cathode-----"
    fig=plt.figure(6)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.75])
    fig.set_size_inches((10., 5.0))
    
    fs = 20
    lw = 2.0
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')    
        
    p1, = plt.plot(SV.loc[:, 'Time'], pct_error_Li, 'k-', linewidth=lw)
    plt.xticks([0, 30000, 60000, 90000, 120000, 150000, 180000])
#    p2, = plt.plot(SV.loc[:, 'Time'], n_Li_cat, 'g-', linewidth=lw)
#    p3, = plt.plot(SV.loc[:, 'Time'], n_Li_Li2S, linewidth=lw)
#    p4, = plt.plot(SV.loc[:, 'Time'], n_Li_tot, linewidth=lw)
#    plt.legend(['Lithium from separator', 'Change in elyte', 'Change in solid', 'Total change'])
    
    fig=plt.figure(7)
    ax = fig.add_axes([0.2, 0.2, 0.6, 0.75])
    fig.set_size_inches((10., 5.0))
    
    fs = 20
    lw = 2.0
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')    
        
    p1, = plt.plot(SV.loc[:, 'Time'], eps_S_vec, 'k-', linewidth=lw)
    p2, = plt.plot(SV.loc[:, 'Time'], eps_Li2S_vec, linewidth=lw)
    p3, = plt.plot(SV.loc[:, 'Time'], eps_el_vec, linewidth=lw)
    p4, = plt.plot(SV.loc[:, 'Time'], eps_C_vec, linewidth=lw)
    plt.legend(['S8', 'Li2S', 'Elyte', 'C'])
    plt.xticks([0, 30000, 60000, 90000, 120000, 150000, 180000])
    
    return

"""========================================================================="""

def plot_sim(tags, SV_df_stage, stage, yax, fig, axes):
    
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
    SV_df = SV_df_stage.copy()
    SV_df.loc[:, 'Time'] *= -cathode.i_ext_amp*inputs.A_cat/3600/(cathode.m_S_0 + cathode.m_S_el)
    print(SV_df.iloc[-1, -1])
    t = SV_df['Time']
    # Plot potential for the electrolyte and the double layer
    SV_plot = SV_df.plot(x='Time', y=phi, ax=axes[0], xlim=[0,t.iloc[-1]])
    SV_plot.set_title(stage, fontsize = fontsize)
    SV_plot.set_ylabel(r'$V_{cell}$ [V]', fontsize = fontsize)
    SV_plot.set_xlabel('Capacity $[A-h/kg_{sulfur}]$', fontsize = fontsize).set_visible(False)
    SV_plot.set_xlim((0, 1750))
    SV_plot.set_ylim((1.5, 2.8))
    SV_plot.set_yticks([1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8])
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
    SV_plot = SV_df.plot(x='Time', y=tags['rho_el'][4:], logy=True, ax=axes[2], xlim=[0,t.iloc[-1]]) #ax=axes[2]
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
#    SV_plot.set_xlim((0, 1750))
#    SV_plot.legend(loc=2, bbox_to_anchor=(1.0, 1), ncol=1, borderaxespad=0,
#                   frameon=False, fontsize = 15).set_visible(False)
#    SV_plot.tick_params(axis='both', labelsize=16)
##    SV_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
#    
#    # Plot species densities in anode electrolyte
#    SV_plot = SV_df.plot(x='Time', y=tags['rho_el_an'][4:], logy=True, ax=axes[4], xlim=[0,t.iloc[-1]]) #
##    SV_plot.set_title(stage, fontsize = fontsize)
#    SV_plot.set_ylabel(r'$\rho_k$ [kmol/m$^3]$', fontsize = fontsize)
#    SV_plot.set_xlabel('Capacity $[Ah/kg_{sulfur}]$', fontsize = fontsize).set_visible(True)
#    SV_plot.set_xlim((0, 1750))
#    SV_plot.legend(loc=2, bbox_to_anchor=(1.0, 1), ncol=1, borderaxespad=0,
#                   frameon=False, fontsize = 15).set_visible(False)
#    SV_plot.tick_params(axis='both', labelsize=16)
##    SV_plot.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    
    
    return

"============================================================================="

def plot_meanPS(SV, tags):
    
    SV_df = SV.copy()
    SV_df.loc[:, 'Time'] *= -cathode.i_ext_amp*inputs.A_cat/3600/(cathode.m_S_0 + cathode.m_S_el)
    
    C_k = SV_df[tags['rho_el'][cathode.i_S8:-2]].copy()
    meanPS = np.zeros([len(C_k.index)])
    for j in np.arange(len(C_k.index)):
        meanPS[j] = sum(cathode.n_S_atoms[0:-2]*C_k.iloc[j, :])/sum(C_k.iloc[j, :])
      
        
        
    "Set up your figure"
    fig = plt.figure(2)
    ax = fig.add_axes([0.2,0.2,0.6,0.75])
    fig.set_size_inches((8.,5.0))
    
    "Formatting for the figure:"
    fs = 20     #font size for plots
    lw = 2.0    #line width for plots
#    font = plt.matplotlib.font_manager.FontProperties(family='Times New Roman',size=fs-1)
    
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fs)
        tick.label1.set_fontname('Times New Roman')    
    
    p1, = plt.plot(SV_df.loc[:, 'Time'], meanPS[:], 'k--', linewidth=lw)
#    p1, = plt.plot(SV_df.loc[:, 'Time'], SV_df.loc[:, tags['phi_ed']], 'k-', linewidth=lw)
    plt.xlim((0, 1770))
    plt.xticks([0, 250, 500, 750, 1000, 1250, 1500, 1750])
#    plt.ylim((1.6, 2.6))
    plt.yticks([2, 3, 4, 5, 6, 7, 8])
    plt.ylabel('Mean PS order', fontstyle='normal', fontname='Times new Roman', fontsize=fs+2, labelpad=5.0)
    plt.xlabel(r'Capacity $[\mathrm{Ah} \hspace{0.5} \mathrm{kg}^{-1}_{\mathrm{sulfur}}]$', fontstyle='normal', fontname='Times new Roman', fontsize=fs+2, labelpad=5.0)
        
    return

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
            newcols_el = {k + offset: 'rho_'+species+'_an'+str(j+1)}
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
            newcols_el = {k + offset: 'rho_'+species+'_sep'+str(j+1)}
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
    

    
    
if __name__ == "__main__":
    conservation_tests(SV_dch, tags)
    
    
    
    
    
    
    