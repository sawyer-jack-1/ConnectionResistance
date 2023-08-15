from utilities.rotation_matrices import Rxyz, Rx2D
from utilities.graphs import ConnectionGraph
from utilities.connection_resistance import ConnectionResistance
import numpy as np
from pathlib import Path
import os

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from definitions import CONFIG_DUMBBELL_GRAPH, CONFIG_1CYCLE_CG, CONFIG_WHEATSTONEBRIDGE_GRAPH
from config.yaml_functions import yaml_loader

if __name__=='__main__':
    #experiment_name = 'dumbell_graph'
    #experiment_name = '1_cycle_graph'
    experiment_name = 'wheatstone_bridge_graph'


    if experiment_name == '1_cycle_graph':
        conf_cg = yaml_loader(CONFIG_1CYCLE_CG)
    elif experiment_name == 'dumbell_graph':
        conf_cg = yaml_loader(CONFIG_DUMBBELL_GRAPH)
    elif experiment_name == 'wheatstone_bridge_graph':
        conf_cg = yaml_loader(CONFIG_WHEATSTONEBRIDGE_GRAPH)

    rotations = np.arange(0, 360, 5)
    rotations[0] = 1
    connectionGraph = ConnectionGraph(conf_cg.copy())
    connectionResistance = ConnectionResistance(conf_cg.copy(), connectionGraph)
    R_standard = connectionResistance.calc_effective_resistance()
    print("R_standard: ", R_standard)

    RC_sync_list = []
    RC_sval_list =[]
    RC_trace_list =[]
    RC_det_list = []
    RC_sub_list = []
    for rotation in rotations:
        edge_rotations = conf_cg['connectionGraph']['edgeRotations']
        #edge_rotations[0][0] = rotation
        #edge_rotations[0][1] = rotation
        edge_rotations[0][0] = rotation
        conf_cg['connectionGraph']['edgeRotations'] = edge_rotations
        connectionGraph = ConnectionGraph(conf_cg)
        connectionResistance = ConnectionResistance(conf_cg, connectionGraph)

        if conf_cg['connectionGraph']['connectionDim'] == 1:
            R_connect_sync = R_standard
        else:
            connectionResistance.optimize()
            R_connect_sync = connectionResistance.calc_connection_resistance_sync()
        Reig, Rtrace, Rdet, Rsub = connectionResistance.calc_connection_resistance_schur()

        RC_sync_list.append(R_connect_sync)
        RC_sval_list.append(Reig)
        RC_trace_list.append(Rtrace)
        RC_det_list.append(Rdet)
        RC_sub_list.append(Rsub)

    RC_sync = np.array(RC_sync_list)
    RC_sval = np.array(RC_sval_list)
    RC_trace = np.array(RC_trace_list)
    RC_det = np.array(RC_det_list)
    RC_sub = np.array(RC_sub_list)
    print("RC_sval: ", RC_sval)

    n = R_standard.shape[0]
    xticks = np.arange(0, 380, 120)
    yticks = np.array([0, 1, 2, 3]) #, 4, 5])
    fig, axes = plt.subplots(n, n, figsize=(25, 10))
    for i in range(0, n):
        for j in range(0, n):
            axes[i, j].plot([0, rotations[-1]], [R_standard[i, j], R_standard[i, j]], label='ER')
            axes[i, j].plot(rotations, RC_sync[:, i, j], label='CR sync')
            #axes[i, j].plot(rotations, RC_sval[:, i, j], linestyle='-.', label='CR sval')
            axes[i, j].plot(rotations, RC_trace[:, i, j], linestyle='--', label='CR trace')
            #axes[i, j].plot(rotations, RC_sub[:, i, j], linestyle='-.', label='CR sub')
            #axes[i, j].plot(rotations, RC_det[:, i, j], linestyle='--', label='CR det')
            axes[i, j].set_xlim([1, 360])
            axes[i, j].set_xticks(xticks)
            #axes[i, j].set_yticks(yticks)
            if i==0 and j==0:
                axes[j, i].legend()
            if j == n-1:
                axes[j, i].set_xlabel('Rotations')
            if i == 0:
                axes[j, i].set_ylabel('Magnitude')
    plt.yscale('log')
    plt.tight_layout()
    folder = os.path.join('Figures', experiment_name)
    Path(folder).mkdir(parents=True, exist_ok=True)
    fig.savefig(os.path.join(folder, 'CR_vs_rot.pdf'), format='pdf')
    fig.savefig(os.path.join(folder, 'CR_vs_rot.png'), format='png')
    plt.show()


