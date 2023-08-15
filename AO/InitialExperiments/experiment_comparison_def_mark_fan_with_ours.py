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
    experiment_name = 'dumbell_graph'
    #experiment_name = '1_cycle_graph'
    #experiment_name = 'wheatstone_bridge_graph'


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

    # Calculate standard effective resistance
    er_standard = connectionResistance.calc_effective_resistance()

    cr_mark_and_fan_def = []
    cr_sync_def=[]
    cr_trace_def = []
    for rotation in rotations:
        edge_rotations = conf_cg['connectionGraph']['edgeRotations']
        edge_rotations[0][0] = rotation

        conf_cg['connectionGraph']['edgeRotations'] = edge_rotations
        connectionGraph = ConnectionGraph(conf_cg)
        connectionResistance = ConnectionResistance(conf_cg, connectionGraph)

        # Calculate syncrhonization connection resistance
        connectionResistance.optimize()
        cr_sync_def.append(connectionResistance.calc_connection_resistance_sync())

        # Calculate trace connection resistance
        cr_trace_def.append(connectionResistance.calc_connection_resistance_trace_def(er_standard))
    cr_mark_and_fan_def = np.array(cr_mark_and_fan_def)
    cr_trace_def = np.array(cr_trace_def)

    n = er_standard.shape[0]
    xticks = np.arange(0, 380, 120)
    yticks = np.array([0, 1, 2, 3])  # , 4, 5])
    fig, axes = plt.subplots(n,n, figsize=(25, 10))
    for i in range(0, n):
        for j in range(0, n):
            axes[i, j].plot([0, rotations[-1]], [er_standard[i, j], er_standard[i, j]], label='Standard ER')
            axes[i, j].plot(rotations, cr_sync_def[:, i, j], label='CR sync')
            axes[i, j].plot(rotations, cr_trace_def[:, i, j], linestyle='--', label='CR trace')
            axes[i, j].set_xlim([1, 360])
            axes[i, j].set_xticks(xticks)
            # axes[i, j].set_yticks(yticks)
            if i == 0 and j == 0:
                axes[j, i].legend()
            if j == n - 1:
                axes[j, i].set_xlabel('Rotations')
            if i == 0:
                axes[j, i].set_ylabel('Magnitude')
    plt.tight_layout()
    folder = os.path.join('Figures', experiment_name)
    Path(folder).mkdir(parents=True, exist_ok=True)
    fig.savefig(os.path.join(folder, 'CR_vs_rot.pdf'), format='pdf')
    fig.savefig(os.path.join(folder, 'CR_vs_rot.png'), format='png')
    plt.show()
