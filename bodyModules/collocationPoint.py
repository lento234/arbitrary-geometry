# -*- coding: utf-8 -*-
"""
Name:           collocationPoint.py
Description:    calcuates the coordinates of collocation points
Author:         Lento Manickathan - 1544101
"""
#==============================================================================
# Importing Modules
#==============================================================================

# Standard scientific module
import numpy as np


#==============================================================================
# Calculates the coordinates of collocation point
#==============================================================================
def collocationPoint(start, end, norm):
    '''
    Calculates the collocation point for given panel coordinates and displaced
    in/out in the 'normal' direction. Collocation Point is used to evalute the 
    boundary conditions at the body: No slip and zero normal velocity.
    
    Input
    -----
    start   
        - x,y panel starting coordinates. 2D array of x,y coordinates in
        row 0 and row 1 respectively.
              
    end
        - x,y panel end coordinates. 2D array of x,y coordinates in row 0 and
        row 1 respectively.
              
    norm
        - x,y unit normal vector component of given x,y coordinates. 2D 
        array of x,y components in row 0 and row 1 respectively. 
              
    Returns
    -------
    collocationPoint 
        - x,y coordinates of collocation point. 2D array of x,y  
        coordinates in row 0 and row 1 respectively.
    '''
    
    # Collocation point
    collocationPoint = (start + end)/2 + norm*np.finfo(np.float64).eps
    
    return collocationPoint
    

#==============================================================================
# Other modules
#==============================================================================
# Calculates the collocation point - [split input]
def collocationPoint_split(x1,x2, y1,y2, xNorm,yNorm):
    '''
    Calculates the collocation point
    '''
    
    xCp = np.add((x1 + x2)/2, xNorm*np.finfo(float).eps) 
    yCp = np.add((y1 + y2)/2, yNorm*np.finfo(float).eps) 
    
    return xCp, yCp