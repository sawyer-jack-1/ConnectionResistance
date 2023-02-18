from utilities.optimizer import ConnectionResistanceOptimizer
from utilities.rotation_matrices import Rxyz, Rx2D, Rx1D
from utilities.schur_comp import schur_comp_er

import numpy as np

def make_init_vector_1D(num_nodes):
    init_vector = []
    for i in range(num_nodes):
        init_vector.append(Rx1D(0))
    return np.array(init_vector)

def make_init_vector_2D(num_nodes):
    init_vector = []
    for i in range(num_nodes):
        init_vector.append(Rx2D(0))
    return np.array(init_vector)

def make_init_vector_3D(num_nodes):
    init_vector = []
    for i in range(num_nodes):
        init_vector.append(Rxyz([0,0,0]))
    return np.array(init_vector)


class ConnectionResistance(ConnectionResistanceOptimizer):
    def __init__(self, config, connection_graph):
        super().__init__(connection_graph)
        self.config = config

        self.connectionDim = self.config['connectionGraph']['connectionDim']
        self.numNodes = self.config['graph']['numNodes']

        if self.connectionDim == 1:
            self.g_init = make_init_vector_1D(self.numNodes)
        if self.connectionDim == 2:
            self.g_init = make_init_vector_2D(self.numNodes)
        if self.connectionDim == 3:
            self.g_init = make_init_vector_3D(self.numNodes)

        self.edges = list(self.connection_graph.get_edges())
        self.connection_laplacian = connection_graph.get_connection_laplacian()
        self.pinv_connection_laplacian = connection_graph.get_pinv_connection_laplacian()
        self.pinv_laplacian = connection_graph.get_pinv_laplacian()

    def optimize(self):
        self.result, self.log = self.solve(self.g_init)
        self.g_opt = self.result.point

    def print_optimization_log(self):
        print("opt log: ", self.log)
        print("g opt: \n")
        print(self.g_opt)

    def form_mij_vector(self, i, j):
        if self.config['optimizer']['g_type'] == 'g_opt':
            g = self.g_opt
        elif self.config['optimizer']['g_type'] == 'g_init':
            g = self.g_init

        istart = i*self.connectionDim
        istop = (i+1)*self.connectionDim
        jstart = j*self.connectionDim
        jstop = (j+1)*self.connectionDim

        m = np.zeros((self.connectionDim*self.numNodes, self.connectionDim))
        m[istart:istop, :] = g[i]
        m[jstart:jstop, :] = - g[j]
        return m

    def calc_connection_resistance_ij(self, i, j):
        mij = self.form_mij_vector(i, j)
        R = np.dot(mij.transpose(), np.dot(self.pinv_connection_laplacian, mij))
        Rij = np.linalg.norm(R, ord=2)
        return Rij

    def calc_connection_resistance_sync(self):
        R = np.zeros((self.numNodes, self.numNodes))
        for i in range(0, self.numNodes):
            for j in range(0, self.numNodes):
                if i != j:
                    R[i,j] = self.calc_connection_resistance_ij(i, j)
        return R

    def calc_schur_comp_ij(self, i, j):
        idx_i = self.connectionDim*i
        idx_j = self.connectionDim*j

        ones = np.ones(self.connectionDim)
        zeros = np.zeros(self.connectionDim)
        Xu = np.diag(np.concatenate((ones, zeros)))[:, 0:self.connectionDim]
        Xv = np.diag(np.concatenate((zeros, ones)))[self.connectionDim::, :].transpose()

        indices = np.array([np.arange(idx_i, idx_i+self.connectionDim), np.arange(idx_j, idx_j+self.connectionDim)]).flatten()
        indices_comp = np.arange(0, self.numNodes*self.connectionDim)
        indices_comp = np.delete(indices_comp, indices, axis=0)
        return schur_comp_er(self.connection_laplacian, Xu, Xv, indices, indices_comp)

    def calc_connection_resistance_schur(self):
        R = np.zeros((self.numNodes, self.numNodes))
        Rsum_sc11sc12 = np.zeros((self.numNodes, self.numNodes))
        Rsub_sc11sc12 = np.zeros((self.numNodes, self.numNodes))
        for i in range(0, self.numNodes):
            for j in range(0, self.numNodes):
                if i != j:
                    R[i, j], Rsub_sc11sc12[i, j] = self.calc_schur_comp_ij(i, j)

        return R, Rsub_sc11sc12

    def calc_effective_resistance(self):
        """Constructs the effective resistance matrix from the psudo inverse of the
        graph laplacian
        :return: Effective resistance matrix as n x n numpy matrix
        """
        n, n = self.pinv_laplacian.shape
        R = np.zeros((n, n))
        for i in range(0, n):
            for j in range(0, n):
                R[i, j] = self.pinv_laplacian[i, i] + self.pinv_laplacian[j, j] - 2 * self.pinv_laplacian[i, j]
        return R


