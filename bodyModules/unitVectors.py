# -*- coding: utf-8 -*-
"""
Name:           unitVectors
Description:    Calculates the unit vectors of the body panels: normal and tangent vector
Author:         Lento Manickathan - 1544101 - lento.manickathan@gmail.com
"""

#==============================================================================
# Importing modules
#==============================================================================

# Standard scientific module
import numpy as np


#==============================================================================
# Calculates the unit NORMAL vector of given panel end and start point
#==============================================================================
def normalVector(start, end):
    ''' 
    Calculates the normal vector of a given panel points
    
    Input
    -----
    start   
        - x,y panel starting coordinates. 2D array of x,y coordinates in
        row 0 and row 1 respectively.
              
    end
        - x,y panel end coordinates. 2D array of x,y coordinates in row 0 and
        row 1 respectively.
              
    Returns
    -------
    norm
        - x,y unit normal vector component of given x,y coordinates. 2D array
        of x,y components in row 0 and row 1 respectively.
    '''
    
    # Length of the panel
    r = ((end[0]-start[0])**2 + (end[1]-start[1])**2)**0.5 #hypotenuse
    
    # Normal vector of the panel
    norm = np.array(((-(end[1]-start[1])/r),
                     ((end[0]-start[0])/r)))  # sinAlpha
    
    return norm

#==============================================================================
# Calculated the unit TANGENT vector of given panel end and start point    
#==============================================================================
def tangentVector(start, end):
    '''
    Calculates the tangent vector of a given panel points
    
    Input
    -----
    start   - x,y panel starting coordinates. 2D array of x,y coordinates in
              row 0 and row 1 respectively.
              
    end     - x,y panel end coordinates. 2D array of x,y coordinates in row
              0 and row 1 respectively.
              
    Returns
    -------
    tang    - x,y unit tangent vector component of given x,y coordinates. 2D 
              array of x,y components in row 0 and row 1 respectively.    
              
    '''
    
    # Length of the panel
    r = ((end[0] - start[0])**2 + (end[1] - start[1])**2)**0.5 #hypothenuse
    
    # Tangent vector of the panel
    tang = np.array((((end[0]-start[0])/r),
                     ((end[1]-start[1])/r))) # cosAlpha, SinAlpha
            
    return tang
    
    
#==============================================================================
# Other modules   
#==============================================================================
def normalVector_split(x1,x2,y1,y2):
    ''' 
    Calculates the normal vector
    '''
        
    # Introducing parameters    
    r = np.sqrt(np.power((x2-x1),2) + np.power((y2-y1),2)) #hypotenuse
    
    normx = -(y2-y1)/r  # sinAlpha
    normy = (x2-x1)/r   # cosAlpha
    
    return normx, normy
 
 # calculate the tangential vector
def tangentVector_split(x1,x2,y1,y2):
    '''
    Calculates the tangent vector
    '''
    
    # Introducing parameters        
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2) #hypothenuse
    
    tangx = (x2 - x1)/r # cosAlpha
    tangy = (y2 - y1)/r # sinAlpha
            
    return tangx, tangy
