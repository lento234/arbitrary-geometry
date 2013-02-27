# -*- coding: utf-8 -*-
"""
Created on Sat Feb 23 14:15:28 2013

@author: Urs
"""

import numpy as np

def vortexInduction(data, collocationPoint=None, also_onSelf = False):       
    '''
    Two-dimensional point vortex:
    Equation (10.9 & 10.10)

    '''
    # number of vortex points
    N = np.shape(data)[1] 
            
    if collocationPoint is None:
        # only evaluted on self
        also_onSelf=True
        M = np.shape(data)[1]
       
        x1 = np.tile(np.array([data[0]]).transpose(), [1,N])
        y1 = np.tile(np.array([data[1]]).transpose(), [1,N])
        
    else:
        # collocation point is given
        
        # also calculate on self
        if also_onSelf is True:
            # Calculate induction on self aswell
            # Total number of points of evaluation is control points and vortex points
            M = np.shape(collocationPoint)[1] + np.shape(data)[1]
            
            x1 = np.tile(np.array([np.concatenate((collocationPoint[0],data[0]))]).transpose(), [1,N])
            y1 = np.tile(np.array([np.concatenate((collocationPoint[1],data[1]))]).transpose(), [1,N])
        
        else:
            # Calculate induction only on control points
            
            # Reshape the data
            M = np.shape(collocationPoint)[1] # number of collocation points
            
            # Collocation Point, evaluation of the vortex on points
            x1 = np.tile(np.array([collocationPoint[0]]).transpose(), [1,N])
            y1 = np.tile(np.array([collocationPoint[1]]).transpose(), [1,N])
        
    # Vortex Coordinates
    x0    = np.tile(data[0], [M,1])
    y0    = np.tile(data[1], [M,1])
    Gamma = np.tile(data[2], [M,1])
    
    
    # Defining preliminary variables
    r = (x1 - x0)**2 + (y1 - y0)**2
    
    # Cut-off radius
    self_induction_points = r==0.
    # removing the data at the location
    r[self_induction_points] = np.NaN
    
    # Equation 10.9 and 10.10
    up = (Gamma/(2*np.pi))*((y1-y0)/r)
    wp = (-Gamma/(2*np.pi))*((x1-x0)/r)
    
    # Updating the cut-off velocity : self-induction is zero
    up[self_induction_points] = 0.    
    wp[self_induction_points] = 0.
    
    
    # Sum to get the total induction
    if collocationPoint is None:
        #u = {'vortex': np.sum(up, axis=1)}
        #w = {'vortex': np.sum(wp, axis=1)}
        
        V_ind = {'vortex': np.array([np.sum(up, axis=1),
                                     np.sum(wp, axis=1)])}
        
    else:
        # there is collocation point
        if also_onSelf is True:
        
            # Total induction
            usum = np.sum(up, axis=1)
            wsum = np.sum(wp, axis=1)

             
            V_ind = {'collocationPoint': np.array([usum[0:np.shape(collocationPoint)[1]],
                                                   wsum[0:np.shape(collocationPoint)[1]]]),
                     'vortex':           np.array([usum[np.shape(collocationPoint)[1]:M],
                                                   wsum[np.shape(collocationPoint)[1]:M]])}
                                                   
        else:
             
            V_ind = {'collocationPoint': np.array([np.sum(up, axis=1),
                                                   np.sum(wp, axis=1)])}
                                                   
    return V_ind
    