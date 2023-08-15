import numpy as np

def nonzero_eigenvalue_product(A, thr=10e-2):
    eval, v = np.linalg.eig(A)

    idx = eval>= thr
    eval_nonzero = eval[eval>= thr]
    wtf = np.prod(eval_nonzero)
    return np.prod(eval_nonzero)