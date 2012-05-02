"""This is the docstring for the calcAdiabat.py module."""

import constants as c
import numpy as np
from new_thermo import findTmoist

def calcAdiabat(press0, thetae0, topPress):
    """
    
    Calculates the temperature-pressure coordinates of a moist adiabat.
    
    Parameters
    - - - - - -
    press0: the initial pressure (Pa)
    thetae0: the equivalent potential temperature (K) of the adiabat
    topPress: the final pressure (Pa)
    
    Returns
    - - - - - -
    coords: array([pressVals, tempVals]): 
                 where pressVals (Pa) and tempVals (K) are 50 x 1 arrays 
                 (i.e. coords is a 50 x 2 array, each row contains a temperature
                  -pressure coordinate)
    
    Tests
    - - - - -
    >>> coords = calcAdiabat(800*100, 300, 1000*100)
    >>> coords.shape
    (50, 2)

    """
    
    pressVals = np.linspace(press0, topPress, 50)
    tempVals = np.zeros(pressVals.shape)
    
    for i in range(pressVals.size):
        tempVals[i] = findTmoist(thetae0, pressVals[i])
    
    return np.column_stack((pressVals, tempVals))

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()

     
         
        
    
    
    
    