import numpy as np

def nudgePress(pressVec):
    """
    
    Returns an array with the same entires as pressVec (1D array), except
    all entries that are equal (within a tolerance) are "nudged" (one of the
    entires is increased by a percentage).
    
    Tests
    - - - -
    >>> p = np.array([1.,1.,2.,3.,3.,4.])
    >>> pnew = nudgePress(p)
    >>> ptest = np.array([1, 1.001, 2, 3, 3.001, 4])
    >>> np.alltrue(abs(ptest - pnew)) < 1.e-8
    True
    
    """
    
    newPress = pressVec
    hit, = np.where(np.abs(np.diff(pressVec)) < 1.e-8)
    newPress[hit+1] = pressVec[hit] + 1.e-3*pressVec[hit]
    return newPress

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
     _test()	