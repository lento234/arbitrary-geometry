# -*- coding: utf-8 -*-
"""
Name:       Bodies. [Body and MultiBody modules]
Description:Contains all the data of body and its function
Author:     Lento Manickathan - 1544101
"""

import numpy as np

# Generates coordinates for a 2D cylinder
def generate2DCylinder(n=100):
    '''
    Generates coordinates for a 2D cylinder (circle)
    
    Input
    -----
    n       - number of discretized panel elements
    
    Returns
    -------
    coordinates - x,y coordinates of the panel points. Closed edge by making 
                  the end point same as the first point. 1D array of x and y 
                  coordinates in row 0 and row 1.
              
    '''
    # Unit Radius, Later multiplied
    R=0.5
    # In polar coordinates
    #r       = np.repeat(R,n+1)  # radius for n+1 point of n panels
    r       = np.repeat(R,n+1) 
    theta   = np.linspace(np.pi,-np.pi,(n+1))
    
    # In cartesian coordinates
    # x0 = 0 and ymid = 0
    coordinates = np.array([r*np.cos(theta)+0.5,
                            r*np.sin(theta)])
    
    
    return coordinates
