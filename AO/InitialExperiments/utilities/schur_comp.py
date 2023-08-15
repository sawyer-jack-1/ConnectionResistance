import numpy as np
from utilities.util import nonzero_eigenvalue_product

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

    positive_det = nonzero_eigenvalue_product(sc)

    Rtrace = 1/np.trace(sc)
    Rdet = 1/(np.power(positive_det, 1 / 3) / 2)
    Rsub = np.linalg.norm(np.linalg.pinv((sc11 - sc12)/2), 2)

    u, s, vh = np.linalg.svd(sc12)
    dim = Xu.shape[1]
    if dim == 1:
        Reig = 1/s[0]
    elif dim > 1:
        Reig = 2/(s[0]+s[1])

    return Reig, Rtrace, Rdet, Rsub
