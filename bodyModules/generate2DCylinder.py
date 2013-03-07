# -*- coding: utf-8 -*-
"""
Name:           generate2DCylinder
Description:    Generates 2D cylinder coordinates for a given number of panels
Author:         Lento Manickathan - 1544101 - lento.manickathan@gmail.com
"""
#==============================================================================
# Importing modules
#==============================================================================

# Standard scientific module
import numpy as np


#==============================================================================
# Generates coordinates for a 2D cylinder 
#==============================================================================
def generate2DCylinder(n=100):
    '''
    Generates coordinates for a 2D cylinder (circle) with a unit diameter of 1.
    
    Input
    -----
    n           
        - number of discretized panel elements. Integer. For n panels n+1
        points are generated; panel starting point and panel end point
        
    
    Returns
    -------
    coordinates 
        - x,y coordinates of the panel points. Closed edge by making the end 
        point same as the first point. 1D array of x and y coordinates in
        row 0 and row 1.
              
    '''
    # Unit Radius, Later multiplied
    R = 0.5
    # In polar coordinates
    #r       = np.repeat(R,n+1)  # radius for n+1 point of n panels
    r       = np.repeat(R,n+1) 
    theta   = np.linspace(np.pi,-np.pi,(n+1))
    
    # In cartesian coordinates
    # x0 = 0 and ymid = 0
    coordinates = np.array([r*np.cos(theta)+0.5,
                            r*np.sin(theta)])
    
    
    return coordinates
