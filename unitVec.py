# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 14:46:59 2012
@description: module contains normVec,tangVec,global2panel, panel2global
@author: lento
"""

from numpy import *
#from scipy import *
#from pylab import *

def normVec(start,end):
    ''' Description: Finds the normal vector of a line
    
    Parameters: starting point,end point 
      array containing x,y coordinates in row 1 and 2 respectively
      
      Shape: 2 by n
      
    Returns: normal vector of line
      array containing normal unit vector,
      x and y componenent in row 1 and 2 respectively
      
      Shape: 2 by n
      
    '''
    
    # Introducing parameters
    x1,y1 = start
    x2,y2 = end
    
    r = sqrt((x2-x1)**2 + (y2-y1)**2) #hypotenuse
    
    sinAlpha = (y2-y1)/r
    cosAlpha = (x2-x1)/r
    
    norm = concatenate(([-sinAlpha],[cosAlpha]))

    return norm

def tangVec(start,end):
    ''' Description: Finds the tangential vector of a line
    
    Parameters: starting point,end point
      array containing x,y coordinates in row 1 and 2 respectively

      Shape: 2 by n

    Returns: tangent vector of line
      array containing tangent unit vector,
      x and y componenent in row 1 and 2 respectively

      Shape: 2 by n
      
    '''
    
    # Introducing parameters
    x1,y1 = start
    x2,y2 = end
    r = sqrt((x2 - x1)**2 + (y2 - y1)**2) #hypothenuse

    cosAlpha = (x2 - x1)/r
    sinAlpha = (y2 - y1)/r

    tang = concatenate(([cosAlpha],[sinAlpha]))
    
    return tang
    
def global2panel(point,origin,end):
    '''Description: Transforms Global coordinates to Panel coordinates, Eq(11.23a)
  
    Parameters: points to transform, assosciating origin point, end point
      array containing x,y coordinates in
    
      Shape: meshgrid (no. of points by no. of panels)
    
    Returns: points in panel coordinate
      array containing x,y coorinates in variable xp and yp
    
      Shape: meshgrid (no. of points by no. of panels)
  
    '''	
    
    # Parameters
    xG, yG = point # point to be transformed
    x1, y1 = origin
    x2, y2 = end
    
    r = sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    cosAlpha = (x2 - x1)/r
    sinAlpha = (y2 - y1)/r
    
    xp = cosAlpha*(xG - x1)  + sinAlpha*(yG - y1)
    yp = -sinAlpha*(xG - x1) + cosAlpha*(yG - y1)
    
    return xp,yp
    
def panel2global(velocityPanel,origin,end):
    '''Description: Transforms panel velocity component to global direction 
    
    Parameters: velocity component in panel axis, origin of panel, end of panel
        array containing x,y component in velocityPanel
        
        Shape: meshgrid (no. of points by no. of panels)
        
    Returns: velocity component in global axis
        array containing u and w globabl velocity component
        
        Shape: meshgrid (no. of points by no. of panels)
    
    '''
 
    # Variable definitions
    up,wp = velocityPanel # can be coordinate or velocity component(if later x=u, y=w)
    x1,y1 = origin
    x2,y2 = end
    
    #Figure 11.17
    r = sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    cosAlpha = (x2 - x1)/r
    sinAlpha = (y2 - y1)/r
    
    #Equation 10.7
    u = cosAlpha*up - sinAlpha*wp
    w = sinAlpha*up + cosAlpha*wp
    
    return u,w
