# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 09:18:39 2019

@author: dkorff
"""

import numpy as np
import cantera as ct
from math import pi
from matplotlib import pyplot as plt
from li_s_battery_inputs import inputs

"Import cantera objects - this step is the same regardless of test type"
elyte = ct.Solution(inputs.ctifile, inputs.elyte_phase)
sulfur = ct.Solution(inputs.ctifile, inputs.cat_phase1)
Li2S = ct.Solution(inputs.ctifile, inputs.cat_phase2)
carbon = ct.Solution(inputs.ctifile, inputs.cat_phase3)
conductor = ct.Solution(inputs.ctifile, inputs.metal_phase)

sulfur_el_s = ct.Interface(inputs.ctifile, inputs.sulfur_elyte_phase,
                             [sulfur, elyte, conductor])
Li2S_el_s = ct.Interface(inputs.ctifile, inputs.Li2S_elyte_phase,
                             [Li2S, elyte, conductor])
carbon_el_s = ct.Interface(inputs.ctifile, inputs.graphite_elyte_phase,
                             [carbon, elyte, conductor])

plt.close('all')
N = 100


elyte.electric_potential = 1.0
V_cell = 2.35
carbon.electric_potential = V_cell
conductor.electric_potential = V_cell
#carbon_el_s.electric_potential = 2.5

C_k_0 = np.array([1.023e1, 
                  1.023e1, 
                  1.024, 
                  1.024, 
                  1.943e-2, 
                  1.821e-4, 
                  3.314e-4, 
                  2.046e-5, 
                  5.348e-10, 
                  8.456e-13])

C_k_mat = np.zeros([len(C_k_0), N])
C_k_mat[0, :] = np.linspace(C_k_0[0], C_k_0[0], N)
C_k_mat[1, :] = np.linspace(C_k_0[1], C_k_0[1], N)
C_k_mat[2, :] = np.linspace(C_k_0[2], C_k_0[2], N)
C_k_mat[3, :] = np.linspace(C_k_0[3], C_k_0[3], N)
C_k_mat[4, :] = np.linspace(C_k_0[4], 0.001*C_k_0[4], N)
C_k_mat[5, :] = np.linspace(0.001*C_k_0[5], 10*C_k_0[5], N)
C_k_mat[6, :] = np.linspace(0.001*C_k_0[6], 10*C_k_0[6], N)
C_k_mat[7, :] = np.linspace(0.001*C_k_0[7], 100*C_k_0[7], N)
C_k_mat[8, :] = np.linspace(C_k_0[8], 1e5*C_k_0[8], N)
C_k_mat[9, :] = np.linspace(C_k_0[9], 1e4*C_k_0[9], N)

X_range = np.linspace(0, 1, N)
X_0 = elyte.X
X = elyte.X
G_elyte = np.zeros([elyte.n_species, N])
dG_L = np.zeros([elyte.n_species, N])
dG_S = np.zeros([elyte.n_species, N])
dG_C = np.zeros([elyte.n_species, N, carbon_el_s.n_reactions])
dG_global = np.zeros([elyte.n_species, N])
k_f_L = np.zeros([elyte.n_species, N])
k_f_S = np.zeros([elyte.n_species, N])
k_f_C = np.zeros([elyte.n_species, N, carbon_el_s.n_reactions])
labels = ['$S_8(l)$', '$S_8^{2-}$', '$S_6^{2-}$', '$S_4^{2-}$', '$S_2^{2-}$', '$S^{2-}$']
styles = ['o', 'v', 's', '*', 'x', 'd']
for j in np.arange(4, elyte.n_species):
    for count in np.arange(0, N):
        C_vec = np.copy(C_k_0)
        C_vec[j] = C_k_mat[j, count]
        elyte.X = C_vec/np.sum(C_vec)
        G_elyte[j, count] = elyte.g
        dG_L[j, count] = Li2S_el_s.delta_gibbs
        dG_S[j, count] = sulfur_el_s.delta_gibbs
        dG_C[j, count, :] = carbon_el_s.delta_gibbs
        dG_global[j, count] = Li2S.g - sulfur.g
        k_f_L[j, count] = Li2S_el_s.net_rates_of_progress
        k_f_S[j, count] = sulfur_el_s.net_rates_of_progress
        k_f_C[j, count, :] = carbon_el_s.net_rates_of_progress
        
#        print(carbon_el_s.forward_rates_of_progress, '\n', carbon_el_s.reverse_rates_of_progress, '\n\n')
    
#    plt.figure(1)
#    plt.plot(C_k_mat[j, :], G_elyte[j, :], label=labels[j-4])
#    plt.legend()
#    plt.title('Electrolyte Gibbs Free Energy \n as a function of species concentration')
#    plt.xlabel(r'$C_k$')
#    ax = plt.gca()
#    ax.set_xscale('log')
    
#    plt.figure(18)
#    plt.plot(C_k_mat[j, :], dG_global[j, :], label=labels[j-4])
#    plt.legend()
#    plt.title(r'$\Delta G_{global}$ as a function of species concentration')
#    plt.xlabel(r'$C_k$')
#    ax = plt.gca()
#    ax.set_xscale('log')
    
#plt.figure(2)
#plt.plot(C_k_mat[-1, :], dG_L[-1, :], label=labels[-1])
#plt.title(r'$\Delta G_{rxn}$ at $Li_2S$-Elyte over $S^{2-}$ concentration')
#
#plt.figure(3)
#plt.plot(C_k_mat[4, :], dG_S[4, :], label=labels[0])
#plt.title(r'$\Delta G_{rxn}$ at S-Elyte over $S_8(l)$ concentration')
#
#plt.figure(4)
#plt.plot(C_k_mat[-1, :], k_f_L[-1, :], label=labels[-1])
#plt.title(r'Net rate of progress at $Li_2S$-Elyte over $S^{2-}$ concentration')
#ax = plt.gca()
#ax.set_yscale('log')
#
#plt.figure(5)
#plt.plot(C_k_mat[4, :], k_f_S[4, :], label=labels[0])
#plt.title(r'Net rate of progress at S-Elyte over $S_8(l)$ concentration')
#ax = plt.gca()
#ax.set_yscale('log')

rxn_labels = ['$S_8(l) <-> S_8^{2-}$', 
              '$S_8^{2-} <-> S_6^{2-}$', 
              '$S_6^{2-} <-> S_4^{2-}$', 
              '$S_4^{2-} <-> S_2^{2-}$', 
              '$S_2^{2-} <-> S^{2-}$']
#for j in np.arange(4, elyte.n_species):
#    for k in np.arange(0, carbon_el_s.n_reactions):
#        plt.figure(j+2)
#        plt.plot(C_k_mat[j, :], dG_C[j, :, k], label=rxn_labels[k])
#        plt.title(r'$\Delta G_{rxn}$ at the carbon interface over ' + labels[j-4])
#        plt.legend()
        
rxn_labels = ['$S_8(s) <-> S_8(l)$',
              '$S_8(l) <-> S_8^{2-}$', 
              '$S_8^{2-} <-> S_6^{2-}$', 
              '$S_6^{2-} <-> S_4^{2-}$', 
              '$S_4^{2-} <-> S_2^{2-}$', 
              '$S_2^{2-} <-> S^{2-}$',
              '$2 Li^{+} + S^{2-} <-> Li_2S$']

q_f = np.zeros([sulfur_el_s.n_reactions+carbon_el_s.n_reactions+Li2S_el_s.n_reactions])
q_r = np.zeros([sulfur_el_s.n_reactions+carbon_el_s.n_reactions+Li2S_el_s.n_reactions])

C_k = np.copy(C_k_0)
#C_k[-1] = 1e-13
elyte.X = C_k/np.sum(C_k)

q_f[0] = sulfur_el_s.forward_rates_of_progress
q_r[0] = sulfur_el_s.reverse_rates_of_progress

q_f[1:carbon_el_s.n_reactions+1] = carbon_el_s.forward_rates_of_progress
q_r[1:carbon_el_s.n_reactions+1] = carbon_el_s.reverse_rates_of_progress

q_f[-1] = Li2S_el_s.forward_rates_of_progress
q_r[-1] = Li2S_el_s.reverse_rates_of_progress

x = np.arange(0, len(q_f))
width = 0.35

fig, ax = plt.subplots()
rects1 = ax.bar(x - width/2, q_f, width, label='forward progress')
rects2 = ax.bar(x + width/2, q_r, width, label='reverse progress')

ax.legend()
ax.set_xticks(x)
ax.set_xticklabels(rxn_labels)
ax.set_yscale('log')
plt.show()

        
#for j in np.arange(4, elyte.n_species):
#    for k in np.array((0, 1, 2, 3, 4)):  #np.arange(0, carbon_el_s.n_reactions):
#        plt.figure(j+8)
#        plt.plot(C_k_mat[j, :], k_f_C[j, :, k], label=rxn_labels[k])
#        plt.title(r'Net rate of progress at the carbon interface over ' + labels[j-4])
#        ax = plt.gca()
#        ax.set_yscale('log')
#        plt.legend()
        
    
plt.show()
    
