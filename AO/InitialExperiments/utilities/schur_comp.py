import numpy as np
from utilities.rotation_matrices import Rxyz

def schur_comp(M, indices, indices_comp):
    A = M[np.ix_(indices, indices)]
    B = M[np.ix_(indices, indices_comp)]
    C = M[np.ix_(indices_comp, indices)]
    D = M[np.ix_(indices_comp, indices_comp)]
    Dinv = np.linalg.inv(D)
    return A - np.dot(B, np.dot(Dinv, C))

def schur_comp_er(M, Xu, Xv, indices, indices_comp):
    sc = schur_comp(M, indices, indices_comp)
    sc12 = np.dot(Xu.transpose(), np.dot(sc, Xv))
    sc11 = np.dot(Xu.transpose(), np.dot(sc, Xu))


    Rsub_sc11sc12 = np.linalg.norm(np.linalg.pinv((sc11 - sc12)/2), 2)

    u, s, vh = np.linalg.svd(sc12)
    dim = Xu.shape[1]
    if dim == 1:
        Rsc = 1/s[0]
    elif dim > 1:
        Rsc = 2/(s[0]+s[1])

    return Rsc, Rsub_sc11sc12
