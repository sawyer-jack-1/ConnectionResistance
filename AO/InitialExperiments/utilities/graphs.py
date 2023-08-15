#import numpy as np
import autograd.numpy as np

from scipy.sparse.linalg import inv
from scipy.sparse import csc_matrix
import networkx as nx
import matplotlib.pyplot as plt
import scipy

from utilities.rotation_matrices import Rxyz, Rx2D, Rx1D

class BaseGraph():
    def __init__(self, config):
        self.config = config.copy()

        self.graph_for_visualization = nx.Graph() # Only for visualization
        self.graph_for_visualization.add_weighted_edges_from(self.config['graph']['edge_info'])

        self.edge_info = self.config['graph']['edge_info']
        self.edges = self.get_edges()
        self.numEdges = len(self.edges)

        self.numNodes = self.config['graph']['numNodes']
        self.nodes = self.get_nodes()

        self.weight_matrix = self.construct_weight_matrix()
        self.degree_matrix = self.construct_degree_matrix(self.weight_matrix)
        self.laplacian = self.construct_laplacian(self.degree_matrix, self.weight_matrix)
        self.pinv_laplacian = self.construct_pinv_laplacian(self.laplacian)

    def get_edges(self):
        edges = []
        for info in self.edge_info:
            edges.append((info[0], info[1]))
        return edges

    def get_nodes(self):
        return np.arange(0, self.numNodes)

    def get_laplacian(self):
        return self.laplacian

    def get_pinv_laplacian(self):
        return self.pinv_laplacian

    def get_edge_weights(self):
        edge_weights = []
        for info in self.edge_info:
            edge_weights.append(info[-1])
        return np.array(edge_weights)

    def construct_weight_matrix(self):
        weight_matrix = np.zeros((self.numNodes, self.numNodes))
        for edge_weight in self.config['graph']['edge_info']:
            i, j, w = edge_weight
            weight_matrix[i, j] = w
            weight_matrix[j, i] = w
        return weight_matrix

    def construct_degree_matrix(self, weight_matrix):
        return np.diag(np.sum(weight_matrix, axis=1))

    def construct_laplacian(self, degree_matrix, weight_matrix):
        return degree_matrix - weight_matrix

    def construct_pinv_laplacian(self, laplacian):
        return scipy.linalg.pinv(laplacian, rcond=10**(-10))

    def draw(self):
        if self.config['graph']['directed']:
            nx.draw_networkx(self.graph_for_visualization, arrows=True, with_labels=True, font_weight='bold')
        else:
            nx.draw_networkx(self.graph_for_visualization, arrows=False, with_labels=True, font_weight='bold')


class ConnectionGraph(BaseGraph):
    def __init__(self, config):
        super().__init__(config)
        self.config = config.copy()
        self.edges = list(self.get_edges())
        self.nodes = list(self.get_nodes())
        self.numEdges = len(self.edges)
        self.numNodes = len(self.nodes)
        self.connectionDim = self.config['connectionGraph']['connectionDim']
        self.edgeRotations = self.config['connectionGraph']['edgeRotations']

        self.connection_incidence_matrix = self.construct_connection_incidence_matrix()
        self.connection_incidence_matrix_markfan = self.construct_connection_incidence_matrix_markfan()
        self.connection_laplacian = self.construct_connection_laplacian()
        self.connection_laplacian_markfan = self.construct_connection_laplacian_markfan()
        self.pinv_connection_laplacian = self.construct_pinv_connection_laplacian()
        self.pinv_connection_laplacian_markfan = self.construct_pinv_connection_laplacian_markfan()
        self.set_edge_labels()

    def get_connectionDim(self):
        return self.connectionDim

    def get_numEdges(self):
        return self.numEdges

    def get_numNodes(self):
        return self.numNodes

    def get_connection_laplacian(self):
        return self.connection_laplacian

    def get_connection_laplacian_markfan(self):
        return self.connection_laplacian_markfan

    def get_connection_incidence_matrix(self):
        return self.connection_incidence_matrix

    def get_pinv_connection_laplacian(self):
        return self.pinv_connection_laplacian

    def get_pinv_connection_laplacian_markfan(self):
        return self.pinv_connection_laplacian_markfan

    def get_connection_incidence_matrix_markfan(self):
        return self.connection_incidence_matrix_markfan

    def get_connection_laplacian_markfan(self):
        return self.connection_laplacian_markfan

    def get_pinv_connection_laplacian_markfan(self):
        return self.pinv_connection_laplacian_markfan

    def get_eigenvectors(self, i=None):
        if i !=None:
            return self.connection_evec[:, i]
        else:
            return self.connection_evec
    
    def eigdecomp_connection_laplacian(self):
        eval, evec = np.linalg.eig(self.connection_laplacian)
        idx = np.argsort(eval)
        #eval_sort = eval[idx]
        #evec_sort = evec[:, idx]
        #wtf = evec[:, 0]
        #norm = np.dot(evec_sort[:, 0].transpose(), evec_sort[:, 0])
        #norm2 = np.dot(evec_sort[:, 10].transpose(), evec_sort[:, 10])
        #norm3 = np.dot(evec[: ,0].transpose(), evec[:, 1])
        return eval[idx], evec[:, idx]

    def set_edge_labels(self):
        edge_labels = {}
        edge_info_graph_for_visualization = nx.get_edge_attributes(self.graph_for_visualization, 'weight')
        edges_graph_for_visualization = self.graph_for_visualization.edges()
        for edge_nr, edge in enumerate(edges_graph_for_visualization):
            weight = edge_info_graph_for_visualization[edge]
            label = self.edgeRotations[edge_nr].copy()
            label.append(weight)
            edge_labels[edge] = label
        nx.set_edge_attributes(self.graph_for_visualization, edge_labels, "edge_labels")

    def construct_connection_incidence_matrix(self):
        connection_incidence_matrix = np.zeros((self.numEdges*self.connectionDim, self.numNodes*self.connectionDim))
        for edge_nr, edge_rot in enumerate(self.edgeRotations):
            e_start = edge_nr * self.connectionDim
            e_stop = (edge_nr+1) * self.connectionDim
            i_start = self.edges[edge_nr][0] * self.connectionDim
            i_stop = (self.edges[edge_nr][0] + 1) * self.connectionDim
            j_start = self.edges[edge_nr][1] * self.connectionDim
            j_stop = (self.edges[edge_nr][1] + 1) * self.connectionDim

            connection_incidence_matrix[e_start:e_stop, i_start:i_stop] = np.diag(np.ones(self.connectionDim))
            if self.connectionDim == 1:
                connection_incidence_matrix[e_start:e_stop, j_start:j_stop] = - Rx1D(self.edgeRotations[edge_nr][0])
            elif self.connectionDim == 2:
                connection_incidence_matrix[e_start:e_stop, j_start:j_stop] = - Rx2D(self.edgeRotations[edge_nr][0])
            elif self.connectionDim == 3:
                connection_incidence_matrix[e_start:e_stop, j_start:j_stop] = - Rxyz(self.edgeRotations[edge_nr])
        return connection_incidence_matrix

    def construct_connection_incidence_matrix_markfan(self):
        connection_incidence_matrix_markfan = np.zeros((self.numEdges * self.connectionDim, self.numNodes * self.connectionDim))
        for edge_nr, edge_rot in enumerate(self.edgeRotations):
            e_start = edge_nr * self.connectionDim
            e_stop = (edge_nr+1) * self.connectionDim
            i_start = self.edges[edge_nr][0] * self.connectionDim
            i_stop = (self.edges[edge_nr][0] + 1) * self.connectionDim
            j_start = self.edges[edge_nr][1] * self.connectionDim
            j_stop = (self.edges[edge_nr][1] + 1) * self.connectionDim

            if self.connectionDim == 1:
                connection_incidence_matrix_markfan[e_start:e_stop, i_start:i_stop] = Rx1D(self.edgeRotations[edge_nr][0])
            elif self.connectionDim == 2:
                connection_incidence_matrix_markfan[e_start:e_stop, i_start:i_stop] = Rx2D(self.edgeRotations[edge_nr][0])
            elif self.connectionDim == 3:
                connection_incidence_matrix_markfan[e_start:e_stop, i_start:i_stop] = Rxyz(self.edgeRotations[edge_nr])
            connection_incidence_matrix_markfan[e_start:e_stop, j_start:j_stop] = - np.diag(np.ones(self.connectionDim))
        return connection_incidence_matrix_markfan

    def construct_connection_laplacian(self):
        edge_weights = self.get_edge_weights()
        edge_weights_extended = np.diag(np.repeat(edge_weights, repeats=self.connectionDim))
        v = np.dot(edge_weights_extended, self.connection_incidence_matrix)
        return np.dot(self.connection_incidence_matrix.transpose(), v)

    def construct_connection_laplacian_markfan(self):
        edge_weights = self.get_edge_weights()
        edge_weights_extended = np.diag(np.repeat(edge_weights, repeats=self.connectionDim))
        v = np.dot(edge_weights_extended, self.connection_incidence_matrix_markfan)
        return np.dot(self.connection_incidence_matrix_markfan.transpose(), v)

    def construct_pinv_connection_laplacian(self):
        return scipy.linalg.pinv(self.connection_laplacian, rcond=10**(-10))

    def construct_pinv_connection_laplacian_markfan(self):
        return scipy.linalg.pinv(self.connection_laplacian_markfan, rcond=10**(-10))

    def get_graph_for_visualization(self):
        return self.graph_for_visualization

    def draw(self, save=False):
        fig, ax = plt.subplots()
        plt.title('Graph layout')
        pos = nx.spring_layout(self.graph_for_visualization)
        nx.draw_networkx(self.graph_for_visualization, pos = pos, arrows=True, with_labels=True,
                         font_weight='bold', ax=ax)
        edge_labels = nx.get_edge_attributes(self.graph_for_visualization, 'edge_labels')
        nx.draw_networkx_edge_labels(self.graph_for_visualization,
                                        font_size = 7,
                                         pos=pos,
                                         edge_labels=edge_labels,
                                         label_pos=0.5,
                                         rotate = True,
                                         ax = ax)
        fig.savefig(f'Figures/graph_with_weight_and_rot.pdf', format='pdf')


from definitions import CONFIG_1CYCLE_CG
from config.yaml_functions import yaml_loader
conf_cg = yaml_loader(CONFIG_1CYCLE_CG)

if __name__ == '__main__':
    baseGraph = BaseGraph(conf_cg)
    connectionGraph = ConnectionGraph(conf_cg)
    connectionGraph.draw()
    plt.show()
