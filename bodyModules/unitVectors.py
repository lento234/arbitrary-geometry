# -*- coding: utf-8 -*-
"""
Name:           Modules for calculating geometry propetries
Description:    Contains all the data of body and its function
Author:         Lento Manickathan - 1544101
"""


import numpy as np

# calculate the normal vectors
def normalVector(start, end):
    ''' 
    Calculates the normal vector
    '''
    
    # Introducing parameters    
    r = np.sqrt(np.power((end[0]-start[0]),2) + np.power((end[1]-start[1]),2)) #hypotenuse
    
    norm = np.array(((-(end[1]-start[1])/r),((end[0]-start[0])/r)))  # sinAlpha
    
    return norm
    
 # calculate the tangential vector
def tangentVector(start, end):
    '''
    Calculates the tangent vector
    '''
    
    # Introducing parameters        
    r = np.sqrt(np.power((end[0] - start[0]),2) + np.power((end[1] - start[1]),2)) #hypothenuse
    
    tang = np.array((((end[0]-start[0])/r),((end[1]-start[1])/r))) # cosAlpha, SinAlpha
            
    return tang
    
# Split Variables
   
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
