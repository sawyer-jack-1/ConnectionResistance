import pymanopt
import pymanopt.manifolds
import pymanopt.optimizers
from pymanopt.manifolds.stiefel import Stiefel
from pymanopt.manifolds.special_orthogonal_group import SpecialOrthogonalGroup
from utilities.rotation_matrices import Rxyz

#import numpy as np
import autograd.numpy as np

import matplotlib.pyplot as plt

from utilities.graphs import ConnectionGraph

class ConnectionResistanceOptimizer():
    def __init__(self, connection_graph):
        self.connection_graph = connection_graph
        self.connectionDim = self.connection_graph.get_connectionDim()
        self.numEdges = self.connection_graph.get_numEdges()
        self.numNodes = self.connection_graph.get_numNodes()
        self.connection_incidence_matrix = connection_graph.get_connection_incidence_matrix()
        self.manifold = self.construct_SO_manifold()

    def construct_SO_manifold(self):
        return SpecialOrthogonalGroup(n=self.connection_graph.get_connectionDim(),
                                      k=self.connection_graph.get_numEdges())

    def solve(self, initial_point=None):

        #cost = lambda x: np.trace(np.dot(np.dot(self.connection_incidence_matrix, x).transpose()), np.dot(self.connection_incidence_matrix, x))
        #grad = lambda x: 2*np.dot(self.connection_incidence_matrix, x)
        #cost = pymanopt.function.autograd(self.manifold)(cost)
        #grad = pymanopt.function.autograd(self.manifold)(grad)

        def dummy_x():
            I = np.diag(np.ones(self.connectionDim))
            x = []
            for i in range(self.numEdges):
                x.append((i+1)*I)
            return np.array(x)

        def debug_grad():
            xdummy = dummy_x()
            x = xdummy.reshape(self.numEdges*self.connectionDim, -1)
            Bx1 = np.dot(self.connection_incidence_matrix, x)
            Bx = np.dot(self.connection_incidence_matrix, x).reshape(self.numEdges, self.connectionDim, -1)

        @pymanopt.function.autograd(self.manifold)
        def grad(x):
            x = x.reshape(self.numEdges*self.connectionDim, -1)
            #debug_grad()
            Bx = np.dot(self.connection_incidence_matrix, x).reshape(self.numEdges, self.connectionDim, -1)
            return 2*Bx

        @pymanopt.function.autograd(self.manifold)
        def cost(x):
            x = x.reshape(self.connectionDim*self.numNodes, -1)
            #Bx = np.dot(self.connection_incidence_matrix, x)
            #t = np.dot(Bx.transpose(), Bx)
            #traceold = np.trace(t)
            eval = np.dot(self.connection_graph.connection_laplacian, x)
            connection_laplacian = self.connection_graph.get_connection_laplacian()
            cost = np.dot(x.transpose(), np.dot(connection_laplacian, x))
            return np.trace(cost)

        #problem = pymanopt.Problem(self.manifold, cost, euclidean_gradient=grad)
        problem = pymanopt.Problem(self.manifold, cost)

        optimizer = pymanopt.optimizers.SteepestDescent(log_verbosity=2)
        result = optimizer.run(problem, initial_point=initial_point)
        log = {k: optimizer._log['iterations'][k] for k in ('iteration', 'cost', 'gradient_norm')}
        return result, log


from definitions import CONFIG_1CYCLE_CG
from config.yaml_functions import yaml_loader
conf_cg = yaml_loader(CONFIG_1CYCLE_CG)

if __name__ == '__main__':
    connectionGraph = ConnectionGraph(conf_cg)
    optimizer = ConnectionResistanceOptimizer(connectionGraph)

    ref1 = Rxyz([0, 0, 0])
    ref2 = Rxyz([90, 0, 0])
    ref3 = Rxyz([180, 0, 0])
    ref4 = Rxyz([270, 0, 0])

    initial_point = np.array([ref1, ref2, ref3, ref4])
    result, log = optimizer.solve(initial_point)
    #result, log = optimizer.solve()

    g = result.point

    print(result.point)

    connectionGraph.draw()

    plt.figure()
    plt.plot(log['iteration'], log['cost'], label='cost')
    plt.legend()

    plt.figure()
    plt.plot(log['iteration'], log['gradient_norm'], label='gradient norm')
    plt.legend()

    plt.figure()
    plt.plot(log['iteration'], log['gradient_norm'], label='gradient norm')
    plt.yscale('log')
    plt.legend()
    plt.show()






