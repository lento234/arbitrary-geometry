# -*- coding: utf-8 -*-
"""
Name:           getCoordinates
Description:    Extracts the coordinates from a given file.
Author:         Lento Manickathan - 1544101 - lento.manickathan@gmail.com
"""
#==============================================================================
# Importing modules
#==============================================================================

# Standard scientific module
import numpy as np


#==============================================================================
# Gets coordinates from a txt file
#==============================================================================
def getCoordinates(filename, clockwise=True):
    '''
    Get coordintes from text file and returns it.
    
    Input
    -----
    filname     
        - dir/filename.extension. Give the file name with directory path
        and filename with extension. Textfile prefered.
                  
    clockwise   
        - Is the geometry described 'clockwise' or anti-clockwise. Clockwise
        is prefered, starting from the trailing-edge.                  
    
    
    Returns
    -------
    coordinates 
        - x and y coordinates from the file. 1D array of coordinates with 
        coordinates described in clockwise fashion. The geometry is
        normalized with the chord length. The geometry is displaced to the
        x0,y0 = 0,0. The leading edge is placed at x=0 and the airfoil is  
        located on the x-axis.
    '''
    
    # Extracting data from the given file
    body = np.fromfile(filename, sep=' ')
    
    # Organizing data into x and y variable
    #   x: starting from x = 0; y: mid-point is on x-axis; body normalized with chord length
    coordinates = np.array([(body[::2]  - np.min(body[::2])),
                            (body[1::2] - (np.max(body[1::2]) + np.min(body[1::2]))/2)])\
                            /np.max(body[::2])
    
    # The geometry needs to be in clockwise order
    if clockwise is False:
        coordinates = coordinates[:,::-1] # Flipping the array        
        
    return coordinates