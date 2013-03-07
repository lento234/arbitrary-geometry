# -*- coding: utf-8 -*-
"""
Name:           bodyTransform.py
Description:    Modules to transform the coordinates from body geometry to global
Author:         Lento Manickathan - 1544101 - lento.manickathan@gmail.com
"""
#==============================================================================
# Importing modules
#==============================================================================

# Standard scientific module
import numpy as np

#==============================================================================
# Transformation coordinates of body system to a higher level (global)
#==============================================================================
def body2global(coordinates, theta, origin=[0.,0.]):
    '''
    Transforms the coordinates from body to global geometry. Rotates the
    coordinates by theta about origin.
    
    Input
    -----
    coordinates     
        - x,y-coordinates of the body. x,y coordinates in in row 0 and row 1 
        respectively. 
    
    theta
        - Degree of rotation of the body around the origin. [deg]
    
    origin
        - origin of the body coordinates. body is translated to zero the 
        origin point. 
                      
    Returns
    -------
    global_coordinates 
        - x,y coordinates of the body in global geometry. Same as 
        'coordinates'. 2D array of x,y coordinates in row 0 and row 1.
    '''
    
    # converting theta into radians
    theta = np.deg2rad(theta)     
    
    # Rotating the coordinates and translating it the origin.    
    global_coordinates = np.array([\
                ( (coordinates[0]*np.cos(theta)) + (coordinates[1]*np.sin(theta)) + origin[0]),
                (-(coordinates[0]*np.sin(theta)) + (coordinates[1]*np.cos(theta)) + origin[1])
                                  ])
             
    return global_coordinates
    
#==============================================================================
# Other transformation functions
#==============================================================================
def body2global_split(x_coordinates, y_coordinates, x_origin,y_origin, theta_deg):
    '''
    Body to global transformation module
    '''
    
    xG = np.add(np.multiply(x_coordinates,np.cos(np.deg2rad(theta_deg))),
                np.multiply(y_coordinates,np.sin(np.deg2rad(theta_deg))))\
                + x_origin
    
    yG = np.add(- np.multiply(x_coordinates,np.sin(np.deg2rad(theta_deg))),
                np.multiply(y_coordinates,np.cos(np.deg2rad(theta_deg))))\
                + y_origin
                
    return xG, yG
    
# transformation module, from global to panel geometry
def global2panel(point, origin, end):
    ''' global2panel(point, origin, end) -> xp, yp
    
    Transforms Global coordinates to Panel coordinates, Eq(11.23a)
  
    Parameters
    ----------
    points to transform, assosciating origin point, end point array containing x,y coordinates in
    
    Shape: meshgrid (no. of points by no. of panels)
    
    Returns
    -------
    points in panel coordinate array containing x,y coorinates in variable xp and yp
    
    Shape: meshgrid (no. of points by no. of panels)  
    '''	
    
    # Parameters
    xG, yG = point # point to be transformed
    x1, y1 = origin
    x2, y2 = end
    
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    cosAlpha = (x2 - x1)/r
    sinAlpha = (y2 - y1)/r
    
    xp = cosAlpha*(xG - x1)  + sinAlpha*(yG - y1)
    yp = -sinAlpha*(xG - x1) + cosAlpha*(yG - y1)
    
    return xp, yp