# -*- coding: utf-8 -*-
"""
Name:           bodyTransform.py
Description:    Main command module. 
Author:         Lento Manickathan - 1544101
"""

import numpy as np

def body2global(coordinates, theta, origin =[0.,0.]):
    '''
    Body to global transformation module
    '''
    
    # converting theta into radians
    theta = np.deg2rad(theta)     
        
    #    global_coordinates = np.array([\
    #    
    #                np.add(np.multiply(x_coordinates,np.cos(theta)), 
    #                       np.multiply(y_coordinates,np.sin(theta)))\
    #                       + x_origin,\
    #                       
    #                np.add(- np.multiply(x_coordinates,np.sin(theta)),
    #                       np.multiply(y_coordinates,np.cos(theta)))\
    #                       + y_origin\
    #                
    #                                ])
    
    global_coordinates = np.array([\
                ( (coordinates[0]*np.cos(theta)) + (coordinates[1]*np.sin(theta)) + origin[0]),
                (-(coordinates[0]*np.sin(theta)) + (coordinates[1]*np.cos(theta)) + origin[1])
                                  ])
             
    return global_coordinates
	
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