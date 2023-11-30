import numpy as np

#======= functions =======#

def transform_positive( x, a ):
    if x <= a:
        return -np.inf
    return np.log( x - a )

def inverse_transform_positive( x, a ):
    return a + np.exp( x )

def transform_bounded( x, a = 0., b = 1. ):
    if ( x < a or x > b ):
        raise ValueError("x is not in range [a,b].")
    if ( x == a ):
        return -np.infty
    elif ( x == b ):
        return np.infty
    else:
        return np.log( ( x + a ) / ( b - x ) ) 

def inverse_transform_bounded( x, a = 0., b = 1. ):
    if x == -np.inf:
        return a
    elif x == np.inf:
        return b
    else:
        return ( a + b * np.exp( x ) ) / ( 1 + np.exp( x ) )

def hard_penalty( x, a = -np.inf, b = np.inf ):
    if ( x < a ) or ( x > b ):
        return -np.infty
    else:
        return 0.