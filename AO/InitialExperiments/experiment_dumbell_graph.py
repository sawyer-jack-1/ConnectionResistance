from utilities.graphs import ConnectionGraph
from utilities.connection_resistance import ConnectionResistance
from utilities.visualization import plot_graph_with_effective_resistance, plot_effective_resistance
import numpy as np

import matplotlib.pyplot as plt
from definitions import CONFIG_DUMBBELL_GRAPH
from config.yaml_functions import yaml_loader
conf_cg = yaml_loader(CONFIG_DUMBBELL_GRAPH)

if __name__ == '__main__':
    experiment_name = 'dumbbell_graph'
    connectionGraph = ConnectionGraph(conf_cg)
    connectionGraph.draw()
    plt.show()

    connectionResistance = ConnectionResistance(conf_cg, connectionGraph)
    connectionResistance.optimize()
    connectionResistance.print_optimization_log()
    R_standard = connectionResistance.calc_effective_resistance()
    R_connect_sync = connectionResistance.calc_connection_resistance_sync()
    R_connect_sval_sc12, R_connect_sub_sc11sc12 = connectionResistance.calc_connection_resistance_schur()

    print("Rconnect: \n", R_connect_sync)
    print("Rstandard: \n", R_standard)
    print("R_schur_comp: \n", R_connect_sval_sc12)

    weights = np.array(conf_cg['graph']['edge_info'])[:, 2]

    difference = R_connect_sync - R_standard
    print("Difference connect sync and standard:\n ", difference)

    difference = R_connect_sval_sc12 - R_standard
    print("Difference connect schur and standard:\n ", difference)

    plot_graph_with_effective_resistance(connectionGraph, R_standard, folder=experiment_name, title=f'ER')
    plot_graph_with_effective_resistance(connectionGraph, R_connect_sync, folder=experiment_name, title=f'CR')
    plot_effective_resistance(R_standard, R_connect_sync, R_connect_sval_sc12, R_connect_sub_sc11sc12, folder=experiment_name)
    plt.show()