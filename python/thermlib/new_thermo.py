"""This is the docstring for the new_thermo.py module. It contains commonly used 
thermo functions/skewT-lnp functions from ATSC 405"""

import numpy as np
import scipy as sp
from constants import constants
import matplotlib.cbook as cbook
from scipy import optimize 

def convertSkewToTemp(xcoord, press, skew):
    """
    convertSkewToTemp(xcoord, press, skew)

    Determines temperature from knowledge of a plotting coordinate
    system and corresponding plot skew.
    
    Parameters
    - - - - - -
    xcoord : int
        X coordinate in temperature plotting coordinates.
    press : float
        Pressure (hPa).
    skew : int
        Skew of a given coordinate system.

    Returns
    - - - -
    Temp : float
        Converted temperature in degC.

    Examples
    - - - - -
    >>> convertSkewToTemp(300, 8.e4, 30)
    638.6934574096806
    
    """
    Temp = xcoord  + skew * np.log(press);
    return Temp

def convertTempToSkew(Temp, press, skew):
    """
    convertTempToSkew(Temp, press, skew)

    Determines the transformed temperature in plotting coordinates.
    
    Parameters
    - - - - - -
    Temp : float
        Temperature (degC)
    press : float
        Pressure (hPa).
    skew : int
        Designated skew factor of temperature.

    Returns
    - - - -
    tempOut : float
        Converted temperature (degC).

    Examples
    - - - - -
    >>> convertTempToSkew(30., 8.e4, 30)
    -308.69345740968055
    
    """
    
    tempOut = Temp - skew * np.log(press);
    return tempOut

def esat(T):
    """
    esat(T)

    Calculates the saturation water vapor pressure over a flat
    surface of water at temperature 'T'.

    Parameters
    - - - - - -
    T : float or array_like
        Temperature of parcel (K).

    Returns
    - - - -
    esatOut : float or list
        Saturation water vapour pressure (Pa).

    Examples
    - - - - -
    >>> esat(300.)
    3534.5196668891358
    >>> es = esat(np.array([300., 310.]))
    >>> estest = [3534.5196668891358, 6235.5321818976754]
    >>> np.all(abs(estest-es) < 1e-8)
    True
    
    References
    - - - - - -
    Emanuel 4.4.14 p. 117
      
    """
    # determine if T has been input as a vector
    is_scalar=True
    if cbook.iterable(T):
        is_scalar = False
    T=np.atleast_1d(T)
    Tc = T - 273.15
    esatOut = 611.2 * np.exp(17.67 * Tc / (Tc + 243.5))
    # if T is a vector
    if is_scalar:
        esatOut = esatOut[0]
    return esatOut

def wsat(Temp, press):
    """
    wsat(Temp, press)

    Calculates the saturation vapor mixing ratio of an air parcel.

    Parameters
    - - - - - -
    Temp : float or array_like
        Temperature in Kelvin.
    press : float or array_like
        Pressure in Pa.

    Returns
    - - - -
    theWs : float or array_like 
        Saturation water vapor mixing ratio in (kg/kg).

    Raises
    - - - -
    IOError
        If both 'Temp' and 'press' are array_like.

    Examples
    - - - - -
    >>> wsat(300, 8e4)
    0.028751159650442507
    >>> wsat([300,310], 8e4)
    [0.028751159650442507, 0.052579529573838296]
    >>> wsat(300, [8e4, 7e4])
    [0.028751159650442507, 0.033076887758679716]
    >>> wsat([300, 310], [8e4, 7e4])
    Traceback (most recent call last):
        ...
    IOError: Can't have two vector inputs.

    """
    c = constants();
    es = esat(Temp);

    # We need to test for all possible cases of (Temp,press)
    # combinations (ie. (vector,vector), (vector,scalar),
    # (scalar,vector), or (scalar,scalar).
    
    try:
        len(es)
    except:
        esIsVect = False
    else:
        esIsVect = True

    try:
        len(press)
    except:
        pressIsVect = False
    else:
        pressIsVect = True
    
    if esIsVect and pressIsVect:
        raise IOError, "Can't have two vector inputs."
    elif esIsVect:
        theWs = [(c.eps * i/ (press - i)) for i in es]
        # Limit ws values so rootfinder doesn't blow up.       
        theWs = list(replaceelem(theWs,0,0.060))
    elif pressIsVect:
        theWs = [(c.eps * es/ (i - es)) for i in press]
        # Limit ws values so rootfinder doesn't blow up.       
        theWs = list(replaceelem(theWs,0,0.060))
    else: # Neither 'es' nor 'press' in a vector.
        theWs = (c.eps * es/ (press - es))
        # Limit ws value so rootfinder doesn't blow up.
        if theWs > 0.060: theWs = 0.060
        elif theWs < 0: theWs = 0

    # Set minimum and maximum acceptable values for theWs.
        
    try:
        len(theWs)
    except:        
        if theWs > 0.060: theWs = 0.060
        elif theWs < 0: theWs = 0
    else:
        theWs = list(replaceelem(theWs, 0, 0.060))

    return theWs


def replaceelem(theList,lowLim,upLim):
    """
    raplaceelem(theList, lowLim, upLim)

    Replaces any elements in 'theList' greater than 'upLim' and less
    than 'lowLim' with the values of 'upLim' and 'lowLim',
    respectively.

    Parameters
    - - - - - -
    theList : array_like
        An array_like structure, the upper and lower bounds of which
        are to be tested.
    lowLim : int
        Number to replace any values within 'theList' that are lower
        than it.
    upLim : int
        Number to replace any values within 'theList' that are higher
        than it.

    Returns
    - - - -
    newList : array
        Augmentation of 'theList' with upper and lower bounds
        accounted for and replaced if necessary.
    """    
    newList = sp.zeros(len(theList))
    for i in theList:
        if i < lowLim:
            newList[theList.index(i)] = lowLim
        elif i > upLim:
            newList[theList.index(i)] = upLim
        else:
            newList[theList.index(i)] = i
    return newList

def theta(*args):
    """
    theta(*args)

    Computes potential temperature.
    Allows for either T,p or T,p,wv as inputs.
    

    Parameters
    - - - - - -
    T : float
        Temperature (K).
    p : float
        Pressure (Pa).


    Returns
    - - - -
    thetaOut : float
        Potential temperature (K).


    Other Parameters
    - - - - - - - - -
    wv : float, optional
        Vapour mixing ratio (kg,kg). Can be appended as an argument
        in order to increase precision of returned 'theta' value.
    
    
    Raises
    - - - -
    NameError
        If an incorrect number of arguments is provided.
    
    
    References
    - - - - - -
    Emanuel p. 111 4.2.11


    Examples
    - - - - -
    >>> theta(300., 8.e4) # Only 'T' and 'p' are input.
    319.72798180767984
    >>> theta(300., 8.e4, 0.001) # 'T', 'p', and 'wv' all input.
    319.72309475657323
    
    """
    c = constants();
    if len(args) == 2:
        wv = 0;
    elif len(args) == 3:
        wv = args[2];
    else:
        raise NameError('need either T,p or T,p,wv');
    
    T = args[0];
    p = args[1];
    power = c.Rd / c.cpd * (1. - 0.24 * wv);
    thetaOut = T * (c.p0 / p) ** power;
    return thetaOut

def thetaep(Td, T, p):
    """
    thetaep(Td, T, p)

    Calculates the pseudo equivalent potential temperature of a
    parcel. 


    Parameters
    - - - - - -
    Td : float
        Dewpoint temperature (K).
    T : float
        Temperature (K).
    p : float
        Pressure (Pa).


    Returns
    - - - -
    thetaepOut : float
        Pseudo equivalent potential temperature (K).


    Notes
    - - -
    Note that the pseudo equivalent potential temperature of an air
    parcel is not a conserved variable.


    References
    - - - - - -
    Emanuel 4.7.9 p. 132


    Examples
    - - - - -
    >>> thetaep(280., 300., 8.e4) # Parcel is unsaturated.
    344.99830738253371
    >>> thetaep(300., 280., 8.e4) # Parcel is saturated.
    321.5302927767795
    
    """
    c = constants();
    if Td < T:
        #parcel is unsaturated
        [Tlcl, plcl] = LCLfind(Td, T, p);
        wv = wsat(Td, p);
    else:
        #parcel is saturated -- prohibit supersaturation with Td > T
        Tlcl = T;
        wv = wsat(T, p);
    
    # $$$   disp('inside theate')
    # $$$   [Td,T,wv]
    thetaval = theta(T, p, wv);
    power = 0.2854 * (1 - 0.28 * wv);
    thetaepOut = thetaval * np.exp(wv * (1 + 0.81 * wv) \
                                   * (3376. / Tlcl - 2.54));
    #
    # peg this at 450 so rootfinder won't blow up
    #
    if(thetaepOut > 450.):
        thetaepOut = 450;
    return thetaepOut

def thetaes(T, p):
    """
    thetaes(T, p)

    Calculates the pseudo equivalent potential temperature of an air
    parcel.

    Parameters
    - - - - - -
    T : float
        Temperature (K).
    p : float
        Pressure (Pa).


    Returns
    - - - -
    thetaep : float
        Pseudo equivalent potential temperature (K).


    Notes
    - - -
    It should be noted that the pseudo equivalent potential
    temperature (thetaep) of an air parcel is not a conserved
    variable.


    References
    - - - - - -
    Emanuel 4.7.9 p. 132


    Examples
    - - - - -
    >>> thetaes(300., 8.e4)
    412.97362667593831
    
    """
    c = constants();
    # The parcel is saturated - prohibit supersaturation with Td > T.
    Tlcl = T;
    wv = wsat(T, p);
    thetaval = theta(T, p, wv);
    power = 0.2854 * (1 - 0.28 * wv);
    thetaep = thetaval * np.exp(wv * (1 + 0.81 * wv) * \
                                (3376. / Tlcl - 2.54))
    #
    # peg this at 450 so rootfinder won't blow up
    #
    if thetaep > 450.:
        thetaep = 450;

    return thetaep

def tinvert_thetae(thetaeVal, wT, p):
    """
    tinvert_thetae(thetaeVal, wT, p)

    Uses a rootfinder to determine the temperature for which the
    pseudo equivilant potential temperature (thetaep) is equal to the
    equivilant potential temperature (thetae) of the parcel.

    Parameters
    - - - - - -
    thetaeVal : float
        Thetae of parcel (K).
    wtotal : float
        Total water mixing ratio (kg/kg).
    p : float
        Pressure of parcel in (Pa).

    Returns
    - - - -
    theTemp : float
        Temperature for which thetaep equals the parcel thetae (K).
    wv : float
        Vapor mixing ratio of the parcel (kg/kg).
    wl : float
        liquid water mixing ratio of the parcel (kg/kg) at 'p'.

    Raises
    - - - -
    IOError
        If 'p' is larger than 100000 Pa.

    Examples
    - - - - -
    >>> tinvert_thetae(300., 0.001, 8.e4)
    
    """
    
    
    if p > 1.e5:
        raise IOError('expecting pressure level less than 100000 Pa')
    # The temperature has to be somewhere between thetae
    # (T at surface) and -40 deg. C (no ice).    
    handle = Tchange
    theTemp = sp.optimize.zeros.brenth(handle, 233.15, \
                                      thetaeVal, (thetaeVal, wT, p));
    [wv,wl] = findWvWl(theTemp, wT, p);
    return theTemp,wv,wl


def Tchange(Tguess, thetaeVal, wT, p):
    [wv, wl] = findWvWl(Tguess, wT, p);
    tdGuess = Tdfind(wv, p);
    # Iterate on Tguess until this function is
    # zero to within tolerance.
    return thetaeVal - thetaep(tdGuess, Tguess, p);

def findTdwv(wv, p):
    """
    
    in:    wv = mixing ratio in kg/kg, p= pressure (Pa)
    
    out:   dewpoint temperature in K
    
    reference: Emanuel 4.4.14 p. 117
    
    """
    c = constants();
    e= wv*p/(c.eps + wv)
    denom=(17.67/log(e/611.2)) - 1.
    Td = 243.5/denom
    Td = Td + 273.15
    
    return Td

def findTmoist(thetaE0, press):
    """
    findTmoist(thetaE0, press)
    
    Calculates the temperatures along a moist adiabat.
    
    Parameters
    - - - - - -
    thetaE0 : float
        Initial equivalent potential temperature (K).
    press : float or array_like
        Pressure (Pa).

    Returns
    - - - -
    Temp : float or array_like
        Temperature (K) of thetaE0 adiabat at 'press'.

    Examples
    - - - - -
    >>> findTmoist(300., 8.e4)
    270.59590841970277
    
    """

    # First determine if press can be indexed
    try: len(press)
    except: #press is a single value
        Temp = optimize.zeros.brenth(thetaEchange, 200, 400, \
                                        (thetaE0, press));
    else: #press is a vector           
        Temp = []
        press = list(press)        
        for i in press:            
            # This assumes that the dewpoint is somewhere between 
            # 250K and 350K.
            Temp.append(optimize.zeros.brenth(thetaEchange, 200, \
                                                 400, (thetaE0, i)));
            #{'in Tmoist: ',i, result(i)}  
        
    return Temp
    

def thetaEchange(Tguess, thetaE0, press):
    """
    thetaEchange(Tguess, thetaE0, press)

    Evaluates the equation and passes it back to brenth.

    Parameters
    - - - - - -
    Tguess : float
        Trial temperature value (K).
    ws0 : float
        Initial saturated mixing ratio (kg/kg).
    press : float
        Pressure (Pa).

    Returns
    - - - -
    theDiff : float
        The difference between the values of 'thetaEguess' and
        'thetaE0'. This difference is then compared to the tolerance
        allowed by brenth.
        
    """
    thetaEguess = thetaes(Tguess, press);
    #{'in change: ',Tguess,press,thetaEguess,thetaE0}
    #when this result is small enough we're done
    theDiff = thetaEguess - thetaE0;
    return theDiff

def Tdfind(wv, p):
    """
    Tdfind(wv, p)

    Calculates the due point temperature of an air parcel.

    Parameters
    - - - - - -
    wv : float
        Mixing ratio (kg/kg).
    p : float
        Pressure (Pa).

    Returns
    - - - -
    Td : float
        Dew point temperature (K).

    Examples
    - - - - -
    >>> Tdfind(0.001, 8.e4)
    253.39429263963504

    References
    - - - - - -
    Emanuel 4.4.14 p. 117
    
    """
    c = constants();    
    e = wv * p / (c.eps + wv);
    denom = (17.67 / np.log(e / 611.2)) - 1.;
    Td = 243.5 / denom;
    Td = Td + 273.15;
    return Td

def LCLfind(Td, T, p):
    """
    LCLfind(Td, T, p)

    Finds the temperature and pressure at the lifting condensation
    level (LCL) of an air parcel.

    Parameters
    - - - - - -
    Td : float
        Dewpoint temperature (K).
    T : float
        Temperature (K).
    p : float
        Pressure (Pa)

    Returns
    - - - -
    Tlcl : float
        Temperature at the LCL (K).
    plcl : float
        Pressure at the LCL (Pa).

    Raises
    - - - -
    NameError
        If the air is saturated at a given Td and T (ie. Td >= T)
    
    Examples
    - - - - -
    >>> [Tlcl, plcl] =  LCLfind(280., 300., 8.e4)
    >>> print [Tlcl, plcl]
    [275.76250387361404, 59518.928699453245]
    >>> LCLfind(300., 280., 8.e4)
    Traceback (most recent call last):
        ...
    NameError: parcel is saturated at this pressure

    References
    - - - - - -
    Emanuel 4.6.24 p. 130 and 4.6.22 p. 129
    
    """
    c = constants();

    hit = Td >= T;
    if hit is True:
        raise NameError('parcel is saturated at this pressure');

    e = esat(Td);
    ehPa = e * 0.01; #Bolton's formula requires hPa.
    # This is is an empircal fit from for LCL temp from Bolton, 1980 MWR.
    Tlcl = (2840. / (3.5 * np.log(T) - np.log(ehPa) - 4.805)) + 55.;

    r = c.eps * e / (p - e);
    #disp(sprintf('r=%0.5g',r'))
    cp = c.cpd + r * c.cpv;
    logplcl = np.log(p) + cp / (c.Rd * (1 + r / c.eps)) * \
              np.log(Tlcl / T);
    plcl = np.exp(logplcl);
    #disp(sprintf('plcl=%0.5g',plcl))

    return Tlcl, plcl


def findWvWl(T, wT, p):
    """
    findWvWl(T, wT, p)

    Computes the vapour and liquid water mixing ratios.

    Parameters
    - - - - - -
    T : float
        Temperature (K).
    wT : float
        Total water mixing ratio (kg/kg).
    p : float
        Pressure (Pa).


    Returns
    - - - -
    wv : float
        Water vapour mixing ratio (kg/kg).
    wl : float
        Liquid water mixing ratio (kg/kg).


    Raises
    - - - -
    AssertionError
        If any of the inputs are in vector form.

    Examples
    - - - - -
    >>> findWvWl(250., 0.01, 8.e4)
    (0.00074331469136857157, 0.0092566853086314283)
    >>> findWvWl(300., 0.01, 8.e4)
    (0.01, 0)
    >>> findWvWl([250.], 0.01, 8.e4)
    Traceback (most recent call last):
        ...
    AssertionError: A vector is not an acceptable input
    
    """
    args = (T, wT, p)
    assert islist(*args) is False , \
           'A vector is not an acceptable input'
    
    wsVal = wsat(T, p)
    if wsVal > wT: #unsaturated
        wv = wT
        wl = 0
    else:  #saturated
        wv = wsVal
        wl = wT - wv
    return wv, wl


def islist(*args):
    """
    Takes any arguments and determines
    if any can be indexed (ie. are lists
    or arrays or tuples). If any are found to be
    indexable, then 'islist' returns 'True'.
    If none are indexable, then 'islist' returns
    'False'.
    """
    truelist = list(np.zeros(len(args)))
    args = list(args)

    count=0

    for i in args:
        try:
            i[0]
        except:
            truelist[count] = False
        else:
            truelist[count] = True
        count += 1
        
    if any(truelist):
        return True
    else:
        return False
    
def nudgePress(pressVec):
    """
    
    Returns an array identical to pressVec (1D array), except all entries that are
    equal (within a tolerance), are "nudged" (one of the entries is 
    increased by a percentage).
    
    Tests
    - - - -
    >>> p = np.array([1.,1.,2.,3.,3.,4.])
    >>> pnew = nudgePress(p)
    >>> ptest = np.array([1, 1.001, 2, 3, 3.001, 4])
    >>> np.alltrue(abs(ptest - pnew)) < 1.e-8
    True
    
    """
    
    newPress = pressVec
    hit = np.nonzero(np.abs(np.diff(pressVec)) < 1.e-8)
    #nonzero returns a tuple containing the indices of the nonzero
    #entries for each dimension
    newPress[hit[0]+1] = pressVec[hit] + 1.e-3*pressVec[hit]
    return newPress
    
    
    

def _test():
    import doctest
    doctest.testmod()

if __name__ == "__main__":
    _test()