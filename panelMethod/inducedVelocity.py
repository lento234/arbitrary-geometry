# -*- coding: utf-8 -*-
"""
Name:       InducedVelocity.py [Body and MultiBody modules]
Description:Contains all the data of body and its function
Author:     Lento Manickathan - 1544101
"""
#==============================================================================
# Importing modules
#==============================================================================
import numpy as np

#==============================================================================
# Calculate the induced velocity due to source
#==============================================================================
def calc(sigma, x0, y0, x1, y1, x2, y2):
    '''
    Calculate the induced velocity due to a source term
    '''
    
    # Transforming from global coordinates to panel coordinates
    
    #Figure 11.17   
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)  

    cosAlpha = (x2 - x1)/r
    sinAlpha = (y2 - y1)/r

    x0mx1 = x0 - x1 # else, calculated twice
    y0my1 = y0 - y1 # else, calculated twice
    
    # In panel coordinates (from global coordinates)    
    x  = cosAlpha*(x0mx1)  + sinAlpha*(y0my1)
    y  = -sinAlpha*(x0mx1) + cosAlpha*(y0my1)    
    
    x2 = cosAlpha*(x2 - x1)+ sinAlpha*(y2 - y1)
     
    # Calculating the induction at panel coordinates
   
    # Defining preliminary variables    
    r1 = x**2 + y**2
    r2 = (x - x2)**2 + y**2
    
    theta1 = np.arctan2(y, x) # theta1
    theta2 = np.arctan2(y, (x - x2))
    
    #Equation 11.21 and 11.22
    up = (sigma/(4*np.pi))*np.log(r1/r2)
    wp = (sigma/(2*np.pi))*(theta2 - theta1)
        
    # Transforming the velocity from panel coordinates to global coordinates
      
    #Equation 10.7
    u = cosAlpha*up - sinAlpha*wp
    w = sinAlpha*up + cosAlpha*wp
    
    return u, w