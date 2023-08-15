from utilities.graphs import ConnectionGraph
from utilities.connection_resistance import ConnectionResistance
import numpy as np
from pathlib import Path
import os

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from utilities.plotstyle_maps import get_color_map, get_dash_map, get_marker_map

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from definitions import CONFIG_MATPLOTSTYLE_PATH
plt.style.use(CONFIG_MATPLOTSTYLE_PATH)

from definitions import CONFIG_WHEATSTONEBRIDGE_GRAPH
from config.yaml_functions import yaml_loader

if __name__=='__main__':
    experiment_name = 'wheatstone_bridge_graph'
    indices = [(0, 1), (1, 2), (2, 3), (3, 0), (1, 3), (0, 2)]
    indices_markfan_exp = [(1, 0), (1, 2), (1, 3), (0, 2)]
    conf_cg = yaml_loader(CONFIG_WHEATSTONEBRIDGE_GRAPH)

    rotations = np.arange(-180, 180, 5)
    rotations = list(rotations)
    rotations.remove(0)
    rotations = np.array(rotations)
    rotations = np.sort(rotations)
    special_rotations = [0]

    connectionGraph = ConnectionGraph(conf_cg.copy())
    connectionResistance = ConnectionResistance(conf_cg.copy(), connectionGraph)

    # Calculate standard effective resistance
    er_standard = connectionResistance.calc_effective_resistance()

    cr_mark_and_fan_def = []
    cr_sync_def = []
    cr_trace_def = []
    cr_offdiag_trace_def = []
    for rotation in rotations:
        edge_rotations = conf_cg['connectionGraph']['edgeRotations']
        edge_rotations[2][0] = rotation

        conf_cg['connectionGraph']['edgeRotations'] = edge_rotations
        connectionGraph = ConnectionGraph(conf_cg)
        connectionResistance = ConnectionResistance(conf_cg, connectionGraph)

        # Calculate syncrhonization connection resistance
        connectionResistance.optimize()
        cr_sync_def.append(connectionResistance.calc_connection_resistance_sync())

        # Calculate mark and fan connection resistance
        cr_mark_fan = connectionResistance.calc_connection_resistance_mark_and_fan()
        cr_mark_and_fan_def.append(connectionResistance.calc_connection_resistance_mark_and_fan())

        # Calculate trace connection resistance
        cr_trace_def.append(connectionResistance.calc_connection_resistance_trace_def(er_standard))
        cr_offdiag_trace_def.append(connectionResistance.calc_connection_resistance_offdiag_trace_def(er_standard))

    cr_mark_and_fan_def_special_rot = []
    for rotation in special_rotations:
        edge_rotations = conf_cg['connectionGraph']['edgeRotations']
        edge_rotations[2][0] = rotation

        conf_cg['connectionGraph']['edgeRotations'] = edge_rotations
        connectionGraph = ConnectionGraph(conf_cg)
        connectionResistance = ConnectionResistance(conf_cg, connectionGraph)

        # Calculate mark and fan connection resistance
        cr_markfan = connectionResistance.calc_connection_resistance_mark_and_fan()
        cr_mark_and_fan_def_special_rot.append(connectionResistance.calc_connection_resistance_mark_and_fan())

    cr_sync_def = np.array(cr_sync_def)
    cr_mark_and_fan_def = np.array(cr_mark_and_fan_def)
    cr_mark_and_fan_def_special_rot = np.array(cr_mark_and_fan_def_special_rot)
    cr_trace_def = np.array(cr_trace_def)
    cr_offdiag_trace_def = np.array(cr_offdiag_trace_def)

    # Formating
    color_map = get_color_map()
    dash_map = get_dash_map()
    marker_map = get_marker_map()

    n = er_standard.shape[0]
    #xticks = np.arange(-180, 185, 60)
    xticks = [-180, -90, 0, 90, 180]
    xtickslabels = [r'$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
    edges = connectionResistance.get_edges()

    # Figure for main
    fig, axes = plt.subplots(2, 3, figsize=(17, 6))
    l=0
    for k, (i,j) in enumerate(indices):
        if k > 2:
            l = 1
        axes[l, k % 3].plot([-180, rotations[-1]], [er_standard[i, j], er_standard[i, j]],
                     color=color_map['color_darkgray'],
                     label='ER')
        axes[l, k % 3].plot(rotations, cr_mark_and_fan_def[:, i, j],
                         dashes=dash_map['dash3'],
                         color=color_map['color_cyan'],
                         label='CR Chung et al (2014)')
        axes[l, k % 3].plot(special_rotations[0], cr_mark_and_fan_def_special_rot[0, i, j],
                         marker='x',
                         markersize='12',
                         color=color_map['color_cyan'])
        axes[l, k % 3].plot(rotations, cr_trace_def[:, i, j],
                     dashes=dash_map['dash4'],
                     color=color_map['color_darkgreen'],
                     label='CR Def. 6.7')
        axes[l, k % 3].set_xlim([-180, 180])
        axes[l, k % 3].set_ylim([0.13, 0.38])
        axes[l, k % 3].set_xticks(xticks)
        axes[l, k % 3].set_xticklabels(xtickslabels)

        if k==0:
            axes[l, k % 3].legend()
            axes[l, k % 3].set_ylabel('Magnitude')
        axes[l, k % 3].set_title(f'Vertices {(i+1, j+1)}')
        axes[l, k % 3].set_xlabel('Rotations')
    plt.tight_layout()
    folder = os.path.join('Figures', experiment_name)
    Path(folder).mkdir(parents=True, exist_ok=True)
    fig.savefig(os.path.join(folder, 'wheatstonebridge_graph_comparisons.pdf'), format='pdf')
    fig.savefig(os.path.join(folder, 'wheatstonebridge_graph_comparisons.png'), format='png')


    # Figure for main text 2
    fig, axes = plt.subplots(1, len(indices_markfan_exp), figsize=(17, 3))
    for k, (i, j) in enumerate(indices_markfan_exp):
        axes[k].plot([0, rotations[-1]], [er_standard[i, j], er_standard[i, j]],
                     color=color_map['color_darkgray'],
                     label='ER')
        axes[k].plot(rotations, cr_mark_and_fan_def[:, i, j],
                         dashes=dash_map['dash3'],
                         color=color_map['color_cyan'],
                         label='CR Chung et al (2014)')
        axes[k].plot(rotations, cr_trace_def[:, i, j],
                     dashes=dash_map['dash3'],
                     color=color_map['color_darkgreen'],
                     label='CR Def 6.7')
        axes[k].set_xlim([0, 360])
        #axes[k].set_ylim([0.5, 1.4])
        axes[k].set_xticks(xticks)

        if k == 0:
            axes[k].legend()
            axes[k].set_ylabel('Effective resistance')
        axes[k].set_title(f'Vertices {(i + 1, j + 1)}')
        axes[k].set_xlabel('Rotation')
    plt.tight_layout()
    folder = os.path.join('Figures', experiment_name)
    Path(folder).mkdir(parents=True, exist_ok=True)
    fig.savefig(os.path.join(folder, 'wheatstonebridge_graph_comparisons_markfan.pdf'), format='pdf')
    fig.savefig(os.path.join(folder, 'wheatstonebridge_graph_comparisons_markfan.png'), format='png')


    # Figure for supplementary
    fig, axes = plt.subplots(1, len(indices), figsize=(17, 3))
    for k, (i,j) in enumerate(indices):
        axes[k].plot([0, rotations[-1]], [er_standard[i, j], er_standard[i, j]],
                     color=color_map['color_darkgray'],
                     label='ER')
        #if (i,j) in edges:
        #    axes[k].plot(rotations, cr_mark_and_fan_def[:, i, j],
        #                dashes=dash_map['dash3'],
        #                color=color_map['color_cyan'],
        #                label='CR (A)')
        axes[k].plot(rotations, cr_trace_def[:, i, j],
                     dashes=dash_map['dash3'],
                     color=color_map['color_darkgreen'],
                     label='CR trace 1')
        axes[k].plot(rotations, cr_offdiag_trace_def[:, i, j],
                     dashes=dash_map['dash4'],
                     color=color_map['color_lightblue'],
                     label='CR trace 2')
        axes[k].plot(rotations, cr_sync_def[:, i, j],
                     dashes=dash_map['dash2'],
                     color=color_map['color_pink'],
                     label='CR sync')
        axes[k].set_xlim([1, 360])
        axes[k].set_ylim([-0.5, 4])
        axes[k].set_xticks(xticks)

        if k==0:
            axes[k].legend()
            axes[k].set_ylabel('Effective resistance')
        axes[k].set_title(f'Vertices {(i+1, j+1)}')
        axes[k].set_xlabel('Rotation')
    plt.tight_layout()
    folder = os.path.join('Figures', experiment_name)
    Path(folder).mkdir(parents=True, exist_ok=True)
    fig.savefig(os.path.join(folder, 'supplement_wheatstone_graph_comparisons.pdf'), format='pdf')
    fig.savefig(os.path.join(folder, 'supplement_wheatstone_graph_comparisons.png'), format='png')

    # Figure grid in supplementary
    yticks = [-3, 0, 3]
    fig, axes = plt.subplots(n,n, figsize=(25, 10), sharex=True, sharey=True)
    for i in range(0, n):
        for j in range(0, n):
            axes[i, j].plot([0, rotations[-1]], [er_standard[i, j], er_standard[i, j]],
                         color=color_map['color_darkgray'],
                         label='ER')
            axes[i, j].plot(rotations, cr_trace_def[:, i, j],
                         dashes=dash_map['dash3'],
                         color=color_map['color_darkgreen'],
                         label='CR trace 1')
            axes[i, j].plot(rotations, cr_offdiag_trace_def[:, i, j],
                         dashes=dash_map['dash4'],
                         color=color_map['color_lightblue'],
                         label='CR trace 2')
            axes[i, j].plot(rotations, cr_sync_def[:, i, j],
                         dashes=dash_map['dash2'],
                         color=color_map['color_pink'],
                         label='CR sync')
            axes[i, j].set_xlim([0, 360])
            axes[i, j].set_xticks([])
            #if i==0 and j ==0:
            #    axes[i, j].legend(loc='center left', bbox_to_anchor=(1, 0.5))
            if i==n-1:
                axes[i, j].set_xticks(xticks)
            if j==0:
                axes[i, j].set_yticks(yticks)

    plt.tight_layout()

    fig.text(0.5, 0.04, 'Rotation', ha='center', fontsize=14)
    fig.text(0.04, 0.5, 'Effective resistance', va='center', rotation='vertical', fontsize=14)
    plt.subplots_adjust(bottom=0.09)
    plt.subplots_adjust(left=0.07)

    folder = os.path.join('Figures', experiment_name)
    Path(folder).mkdir(parents=True, exist_ok=True)
    fig.savefig(os.path.join(folder, 'supplement_grid_wheatstone_graph_comparisons.pdf'), format='pdf')
    fig.savefig(os.path.join(folder, 'supplement_grid_wheatstone_graph_comparisons.png'), format='png')
    plt.show()
