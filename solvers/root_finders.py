## Own code for root finding
## compatible with NUMBA

import numpy as np
from numba import jit

@jit
def brent(f, x1, x3):
    '''
    Root finding from interval without derivatives
    Based on https://mathworld.wolfram.com/BrentsMethod.html
    '''
    assert f(x1) * f(x3) < 0., "Limit points must have different signs."
    
    x2 = (x1 + x3) / 2
    
    niter = 0
    while np.abs(f(x2)) > 1e-3:
        R = f(x2) / f(x3)
        S = f(x2) / f(x1)
        T = f(x1) / f(x3)
    
        P = S * (T*(R-T)*(x3-x2) - (1-R)*(x2-x1))
        Q = (T-1) * (R-1) * (S-1)
        x2 = x2 + P/Q
        
        niter = niter + 1
        if niter > 1000:
            break

    return x2

@jit
def newton(f, x0, y=0, args=None):
    '''
    Root solving using Newton's method
    Can specify f(x)=y rather than f(x)=0
    '''

    a = 1e-3
    
    if args != None:
        for i in range(1000):
            df = (f(x0+a, *args) - f(x0-a, *args)) /2/a
            x0 = x0 - (f(x0, *args)-y) / df
    else:
        for i in range(1000):
            df = (f(x0+a) - f(x0-a)) /2/a
            x0 = x0 - (f(x0)-y) / df
    
    return x0
