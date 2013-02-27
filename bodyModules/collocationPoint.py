# -*- coding: utf-8 -*-
"""
Name:           collocationPoint.py
Description:    Main command module. 
Author:         Lento Manickathan - 1544101
"""
#==============================================================================
# Importing Modules
#==============================================================================

# Standard scientific module
import numpy as np

#==============================================================================
# Functions
#==============================================================================
# Calculates the collocation point
def collocationPoint(start, end, norm):
    '''
    Calculates the collocation point
    '''
    #start, end, norm
    
    cp = np.add((start + end)/2, 100*norm*np.finfo(np.float64).eps) 
    
    return cp
    
# Calculates the collocation point - [split input]
def collocationPoint_split(x1,x2, y1,y2, xNorm,yNorm):
    '''
    Calculates the collocation point
    '''
    
    xCp = np.add((x1 + x2)/2, 100*xNorm*np.finfo(float).eps) 
    yCp = np.add((y1 + y2)/2, 100*yNorm*np.finfo(float).eps) 
    
    return xCp, yCp