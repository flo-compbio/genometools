import numpy as np

def _reorder(E,A):
    assert E.shape[1] == A.shape[1]
    n = E.shape[1]
    for j in range(n):
        E[:,j] = E[A[:,j],j]
    return E

def quantile_normalize(E):
    E = E.copy()
    A = np.argsort(E,axis=0)
    A_inv = np.argsort(A,axis=0)
    
    E = _reorder(E,A)
    mean = np.mean(E,axis=1)
    n = E.shape[1]
    E = np.tile(mean,(n,1)).T
    E = _reorder(E,A_inv)
    return E
