import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
import os
import pandas as pd
from pathlib import Path

def plot_graph_with_effective_resistance(graph, effective_resistance, folder, title):
    graph_for_visualization = graph.get_graph_for_visualization()

    pos = nx.spring_layout(graph_for_visualization)
    edges =graph.get_edges()
    edge_labels = {}
    for edge in edges:
        edge_labels[edge] = np.round(effective_resistance[edge[0], edge[1]], 5)

    fig, ax = plt.subplots()
    plt.title(title)
    nx.draw_networkx(graph_for_visualization, pos=pos, arrows=True, with_labels=True,
                     font_weight='bold', ax=ax)
    nx.draw_networkx_edge_labels(graph_for_visualization,
                                 font_size=7,
                                 pos=pos,
                                 edge_labels=edge_labels,
                                 label_pos=0.5,
                                 rotate=True,
                                 ax=ax)
    folder = os.path.join('Figures', folder)
    Path(folder).mkdir(parents=True, exist_ok=True)
    fig.savefig(os.path.join(folder, f'{title}.pdf'), format='pdf')
    fig.savefig(os.path.join(folder, f'{title}.png'), format='png')

def plot_effective_resistance(standardER, Rsync, Reig, Rtrace, Rdet, Rsub, folder):
    num_nodes = standardER.shape[0]
    fig, ax = plt.subplots(2, 4, figsize=(15, 7))
    i=0
    j=0
    for node_nr in range(0, num_nodes):
        if i > 0 and i % 4 == 0:
            j += 1
            i = 0
        df = pd.DataFrame({'ER': standardER[node_nr, :],
                           'CRsync': Rsync[node_nr, :],
                          'Reig': Reig[node_nr, :],
                          'Rtrace': Rtrace[node_nr, :],
                          #'Rdet': Rdet[node_nr, :],
                          'Rsub': Rsub[node_nr, :]
                           }
                          )

        df.plot.bar(rot=0, title=f'Node {node_nr}', ax=ax[j, i])

        if j == 1:
            ax[j, i].set_xlabel('Node nr')
        if i == 0:
            ax[j, i].set_ylabel('Magnitude')
        i += 1

    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))

    folder = os.path.join('Figures', folder)
    Path(folder).mkdir(parents=True, exist_ok=True)
    fig.savefig(os.path.join(folder, 'CR_vs_ER.pdf'), format='pdf')
    fig.savefig(os.path.join(folder, 'CR_vs_ER.png'), format='png')
