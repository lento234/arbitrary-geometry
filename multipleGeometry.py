# -*- coding: utf-8 -*-
"""
Case:Arbitrary Geometry\n
Methodology: Solving any geometry using Panel Method\n
Description: Contains the arbitrary.py modules\n 
@author: lento
"""

# external modules

from numpy import *
#from scipy import *

############################################################################

# module to calculate the source term (only)
def sourceTerm(geometries):
    
    # Initially arranging and organizing the data
    data = organizeData(geometries)
    
    controlPoints, panelEnd, panelStart, unitNorm, unitTang = reshapeData(data)
    
    print controlPoints.x, controlPoints.y
    print panelEnd.x, panelEnd.y
    print panelStart.x, panelStart.y
    print unitNorm.x, unitNorm.y
    print unitTang.x, unitTang.y
    
    return data

############################################################################
# List of classes 
class dataClass(object):
    '''Generates an empty class    
    '''
    def __init__(self):
        pass

class Geometry(object):
    '''Geometry Class. Defines the name, x and y coordinates
    
    Parameters
    ----------
    Geometry(object) -> 
        object.x, object.y
    '''
    def __init__(self,x,y):
        self.x = x
        self.y = y

# List of functions
def calc_unitVectors(object):
    '''Calculates the unit vectors of the geometry (object)
    
    Parameters
    ----------    
    calc_unitVectors(object) -> 
        object.unitVectors.norm, object.unitVectors.tang
    '''
    
    norm = normVec((object.points[0,:-1],object.points[1,:-1]),(object.points[0,1:],object.points[1,1:]))
    tang = tangVec((object.points[0,:-1],object.points[1,:-1]),(object.points[0,1:],object.points[1,1:]))
    
    object.unitVectors.norm = norm 
    object.unitVectors.tang = tang
    
    return object
    
def calc_controlPoints(object):
    '''Generates the control points
    
    Parameters
    ----------    
    calc_controlPoints(object) ->
        object.controlPoints
    '''
   
    #cpx = (object.points[0][1:] + object.points[0][:-1])/2 + 100*object.unitVectors.norm[0]*finfo(float).eps
    #cpy = (object.points[1][1:] + object.points[1][:-1])/2 + 100*object.unitVectors.norm[1]*finfo(float).eps
    
    cp = add((object.points[:,1:] + object.points[:,:-1])/2, 100*object.unitVectors.norm*finfo(float).eps) # control points pushed out by delta eps
    object.controlPoints = cp
    
    return object
    
def calc_panelPoints(object):
    '''Generates the panel starting and end points
    
    Parameters
    ----------    
    calc_panelPoints(object) -> 
        object.panel.start, object.panel.end
    '''
    
    object.panel.start = object.points[:,:-1]
    object.panel.end = object.points[:,1:]
    
    return object

def organizeData(object):
    '''Method to organize the inital data\n
    
    Parameters
    ----  
    organizeData(geometries) -> 
        data (.points, .unitVectors, .controlPoints, .panel (.start, .end))
    '''
    data = {}
    for name in object:
        data[name] = dataClass() # empty data class
        
        data[name].__setattr__('points',object[name]) # saving x and y data points into .points  
        
        # empty classes
        data[name].__setattr__('unitVectors',dataClass())
        data[name].__setattr__('controlPoints',dataClass())
        data[name].__setattr__('panel',dataClass())
        
        calc_unitVectors(data[name]) # calculate unitVectors
        calc_controlPoints(data[name]) # calculate control points
        calc_panelPoints(data[name])
        
    return data
    
def reshapeData(object):
    '''Reshape the data and concatenate to one matrix, tile to right format
    
    Parameters
    ----------
    reshapeData(data) ->
        controlPoints, panelEnd, panelStart, unitNorm, unitTang
    
    '''
    
    controlPoints   = Geometry(array([]),array([]))
    panelStart      = Geometry(array([]),array([]))
    panelEnd        = Geometry(array([]),array([]))
    unitNorm        = Geometry(array([]),array([]))
    unitTang        = Geometry(array([]),array([]))
    
    # concatenating all the data to one matrix
    for num in range(len(object)):
        controlPoints.x = concatenate((controlPoints.x, object[object.keys()[num]].controlPoints[0]),axis=0)
        controlPoints.y = concatenate((controlPoints.y, object[object.keys()[num]].controlPoints[1]),axis=0)
    
        panelStart.x = concatenate((panelStart.x, object[object.keys()[num]].panel.start[0]),axis=0)    
        panelStart.y = concatenate((panelStart.y, object[object.keys()[num]].panel.start[1]),axis=0)

        panelEnd.x = concatenate((panelEnd.x, object[object.keys()[num]].panel.end[0]),axis=0)    
        panelEnd.y = concatenate((panelEnd.y, object[object.keys()[num]].panel.end[1]),axis=0)
    
        unitNorm.x = concatenate((unitNorm.x, object[object.keys()[num]].unitVectors.norm[0]),axis=0)
        unitNorm.y = concatenate((unitNorm.y, object[object.keys()[num]].unitVectors.norm[1]),axis=0)
    
        unitTang.x = concatenate((unitTang.x, object[object.keys()[num]].unitVectors.tang[0]),axis=0)
        unitTang.y = concatenate((unitTang.y, object[object.keys()[num]].unitVectors.tang[1]),axis=0)
        
    N = len(controlPoints.x) # total number of control points
    M = len(panelStart.x) # total number of panel points
    
    controlPoints.x = tile(array([controlPoints.x]).transpose(),[1,M])
    controlPoints.y = tile(array([controlPoints.y]).transpose(),[1,M])
    
    panelStart.x = tile(panelStart.x,[N,1])
    panelStart.y = tile(panelStart.y,[N,1])
    
    panelEnd.x = tile(panelEnd.x,[N,1])
    panelEnd.y = tile(panelEnd.y,[N,1])
    
    unitNorm.x = tile(array([unitNorm.x]).transpose(),[1,M])
    unitNorm.y = tile(array([unitNorm.y]).transpose(),[1,M])
    
    unitTang.x = tile(array([unitTang.x]).transpose(),[1,M])
    unitTang.y = tile(array([unitTang.y]).transpose(),[1,M])
       
    return controlPoints, panelEnd, panelStart, unitNorm, unitTang

        




############################################################################
# module to calculate the source term
def sor2D(sigma,controlPoint,panelStart,panelEnd):
    '''Description: Constant Strength Source Method\n
    
    Parameters: source strength, control points, panel starting point, end point
        array of x,y coordinates in row 1 and 2 respectively
        source strength in row 1, length of panel numbers
    
    Returns: u and w velocity in meshgrid
    
    Shape n by m
    
    '''
               
    # In panel coordinates (from global coordinates)
    x,y     = global2panel(controlPoint,panelStart,panelEnd) # control point
    x1,y1   = global2panel(panelStart,panelStart,panelEnd) # panel origin 
    x2,y2   = global2panel(panelEnd,panelStart,panelEnd) # panel end point 
    
    # Defining preliminary variables    
    r1  = (x - x1)**2 + (y)**2 # r1²
    r2  = (x - x2)**2 + (y)**2 # r2²
    
    theta1  = arctan2(y, (x - x1)) # theta1
    theta2  = arctan2(y, (x - x2)) # theta2
    
    #Equation 11.21 and 11.22
    up  = (sigma/(4*pi))*log(r1/r2)
    wp  = (sigma/(2*pi))*(theta2 - theta1)
    
    u,w = panel2global((up,wp),panelStart,panelEnd)
       
    return u,w
    
# calculate the normal vector
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

# calculate the tangential vector
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

# transformation module, from global to panel geometry
def global2panel(point,origin,end):
    '''Transforms Global coordinates to Panel coordinates, Eq(11.23a)
  
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
    
    r = sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    cosAlpha = (x2 - x1)/r
    sinAlpha = (y2 - y1)/r
    
    xp = cosAlpha*(xG - x1)  + sinAlpha*(yG - y1)
    yp = -sinAlpha*(xG - x1) + cosAlpha*(yG - y1)
    
    return xp,yp
    
# to transfer velocity from panel to global geometry
def panel2global(velocityPanel,origin,end):
    '''Transforms panel velocity component to global direction 
    
    Parameters
    ----------
    velocity component in panel axis, origin of panel, end of panel array containing x,y component in velocityPanel
        
    Shape: meshgrid (no. of points by no. of panels)
        
    Returns
    -------
    velocity component in global axis array containing u and w globabl velocity component
        
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
