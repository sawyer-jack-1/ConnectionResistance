from utilities.optimizer import ConnectionResistanceOptimizer
from utilities.rotation_matrices import Rxyz, Rx2D, Rx1D
from utilities.schur_comp import schur_comp_er, schur_comp

import numpy as np
from scipy import linalg

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
        self.connection_incidence_matrix = connection_graph.get_connection_incidence_matrix()
        self.connection_incidence_matrix_markfan = connection_graph.get_connection_incidence_matrix_markfan()
        self.connection_laplacian = connection_graph.get_connection_laplacian()
        self.connection_laplacian_markfan = connection_graph.get_connection_laplacian_markfan()
        self.pinv_connection_laplacian = connection_graph.get_pinv_connection_laplacian()
        self.pinv_connection_laplacian_markfan = connection_graph.get_pinv_connection_laplacian_markfan()
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

    def calc_connection_resistance_sync_ij(self, i, j):
        mij = self.form_mij_vector(i, j)
        R = np.dot(mij.transpose(), np.dot(self.pinv_connection_laplacian, mij))
        Rij = np.linalg.norm(R, ord=2)
        return Rij

    def calc_connection_resistance_sync(self):
        R = np.zeros((self.numNodes, self.numNodes))
        for i in range(0, self.numNodes):
            for j in range(0, self.numNodes):
                if i != j:
                    R[i,j] = self.calc_connection_resistance_sync_ij(i, j)
        return R

    def calc_connection_resistance_mark_and_fan_ij(self, psi, i, j):
        psi_ee = psi[i*self.connectionDim:(i+1)*self.connectionDim, j*self.connectionDim:(j+1)*self.connectionDim]
        r_ij = np.linalg.norm(psi_ee, 2)
        return r_ij

    def calc_connection_resistance_mark_and_fan(self):
        wtf1 = self.connection_incidence_matrix
        wtf2 = self.connection_incidence_matrix_markfan
        psi = np.dot(np.dot(self.connection_incidence_matrix_markfan, self.pinv_connection_laplacian_markfan),
                     self.connection_incidence_matrix_markfan.transpose())
        edges = self.connection_graph.get_edges()

        R = np.zeros((self.numNodes, self.numNodes))
        for k in range(0, len(edges)):
            i, j = edges[k]
            psi_ee = psi[k*self.connectionDim:(k+1)*self.connectionDim, k*self.connectionDim:(k+1)*self.connectionDim]
            R[i, j] = np.linalg.norm(psi_ee, 2)
            R[j, i] = np.linalg.norm(psi_ee, 2)

        #for i in range(0, self.numNodes):
        #    for j in range(0, self.numNodes):
        #        if (i,j) in self.edges or (j,i) in self.edges:
        #            R[i,j] = self.calc_connection_resistance_mark_and_fan_ij(psi, i, j)
        return R

    def calc_connection_resistance_trace_def_ij(self, er_standard, i, j):
        idx_i = self.connectionDim * i
        idx_j = self.connectionDim * j

        ones = np.ones(self.connectionDim)
        zeros = np.zeros(self.connectionDim)
        Xu = np.diag(np.concatenate((ones, zeros)))[:, 0:self.connectionDim]
        Xv = np.diag(np.concatenate((zeros, ones)))[self.connectionDim::, :].transpose()

        indices = np.array(
            [np.arange(idx_i, idx_i + self.connectionDim), np.arange(idx_j, idx_j + self.connectionDim)]).flatten()
        indices_comp = np.arange(0, self.numNodes * self.connectionDim)
        indices_comp = np.delete(indices_comp, indices, axis=0)

        sc = schur_comp(self.connection_laplacian, indices, indices_comp)
        c_ii = np.dot(Xu.transpose(), np.dot(sc, Xu))
        c_jj = np.dot(Xv.transpose(), np.dot(sc, Xv))

        c_ii_inv = np.linalg.inv(c_ii)
        c_jj_inv = np.linalg.inv(c_jj)

        basis = linalg.null_space(self.connection_laplacian)
        n, rho = basis.shape # rho is dimension of nullspace
        #r_ij = rho/self.connectionDim*er_standard + 1/(2*self.connectionDim)*np.trace(c_ii_inv + c_jj_inv)

        #r_ij = er_standard + 1/(2*self.connectionDim)*np.trace(c_ii_inv + c_jj_inv)
        r_ij = 1 / (2 * self.connectionDim) * np.trace(c_ii_inv + c_jj_inv)
        return r_ij

    def calc_connection_resistance_trace_def(self, er_standard):
        cr = np.zeros((self.numNodes, self.numNodes))
        for i in range(0, self.numNodes):
            for j in range(0, self.numNodes):
                if i != j:
                    er_standard_ij = er_standard[i, j]
                    cr[i, j] = self.calc_connection_resistance_trace_def_ij(er_standard_ij, i, j)
        return cr


    def calc_connection_resistance_offdiag_trace_def_ij(self, er_standard, i, j):
        idx_i = self.connectionDim * i
        idx_j = self.connectionDim * j

        ones = np.ones(self.connectionDim)
        zeros = np.zeros(self.connectionDim)
        Xu = np.diag(np.concatenate((ones, zeros)))[:, 0:self.connectionDim]
        Xv = np.diag(np.concatenate((zeros, ones)))[self.connectionDim::, :].transpose()

        indices = np.array(
            [np.arange(idx_i, idx_i + self.connectionDim), np.arange(idx_j, idx_j + self.connectionDim)]).flatten()
        indices_comp = np.arange(0, self.numNodes * self.connectionDim)
        indices_comp = np.delete(indices_comp, indices, axis=0)

        sc = schur_comp(self.connection_laplacian, indices, indices_comp)
        c_ij = np.dot(Xu.transpose(), np.dot(sc, Xv))
        c_jj = np.dot(Xv.transpose(), np.dot(sc, Xv))

        c_ij_inv = np.linalg.inv(c_ij)
        c_jj_inv = np.linalg.inv(c_jj)

        r_ij = er_standard + 1/(2*self.connectionDim)*np.trace(c_ij_inv + c_jj_inv)
        return r_ij

    def calc_connection_resistance_offdiag_trace_def(self, er_standard):
        cr = np.zeros((self.numNodes, self.numNodes))
        for i in range(0, self.numNodes):
            for j in range(0, self.numNodes):
                if i != j:
                    er_standard_ij = er_standard[i, j]
                    cr[i, j] = self.calc_connection_resistance_offdiag_trace_def_ij(er_standard_ij, i, j)
        return cr


    def calc_connection_resistance_singval_def_ij(self, i, j):
        idx_i = self.connectionDim * i
        idx_j = self.connectionDim * j

        ones = np.ones(self.connectionDim)
        zeros = np.zeros(self.connectionDim)
        Xu = np.diag(np.concatenate((ones, zeros)))[:, 0:self.connectionDim]
        Xv = np.diag(np.concatenate((zeros, ones)))[self.connectionDim::, :].transpose()

        indices = np.array(
            [np.arange(idx_i, idx_i + self.connectionDim), np.arange(idx_j, idx_j + self.connectionDim)]).flatten()
        indices_comp = np.arange(0, self.numNodes * self.connectionDim)
        indices_comp = np.delete(indices_comp, indices, axis=0)

        c = schur_comp(self.connection_laplacian, indices, indices_comp)
        c_ij = np.dot(Xu.transpose(), np.dot(c, Xv))
        u, s, vh = np.linalg.svd(c_ij)
        r_ij = 2/(s[0]+s[1])
        return r_ij

    def calc_connection_resistance_singval_def(self):
        cr = np.zeros((self.numNodes, self.numNodes))
        for i in range(0, self.numNodes):
            for j in range(0, self.numNodes):
                if i != j:
                    cr[i, j] = self.calc_connection_resistance_singval_def_ij(i, j)
        return cr

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
        Reig = np.zeros((self.numNodes, self.numNodes))
        Rtrace = np.zeros((self.numNodes, self.numNodes))
        Rdet = np.zeros((self.numNodes, self.numNodes))
        Rsub = np.zeros((self.numNodes, self.numNodes))
        for i in range(0, self.numNodes):
            for j in range(0, self.numNodes):
                if i != j:
                    Reig[i, j], Rtrace[i, j], Rdet[i, j], Rsub[i, j] = self.calc_schur_comp_ij(i, j)

        return Reig, Rtrace, Rdet, Rsub

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

    def get_edges(self):
        return self.edges
