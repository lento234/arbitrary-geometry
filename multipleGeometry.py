# -*- coding: utf-8 -*-
"""
Case:Arbitrary Geometry\n
Methodology: Solving any geometry using Panel Method\n
Description: Contains the arbitrary.py modules\n 
@author: lento
"""

# extra modules

from numpy import *
#from scipy import *

############################################################################

def sourceTerm(geometries):
    
    data = {}
    for name in geometries:
        data[name] = dataClass() # empty data class
        # saving x and y data points into .points                
        data[name].__setattr__('points',Geometry('Geometry points',geometries[name][0],geometries[name][1]))
        data[name].__setattr__('unitVectors',dataClass())
        data[name].__setattr__('controlPoints',dataClass())
        data[name].__setattr__('panel',dataClass())
        
        calc_unitVectors(data[name]) # calculate unitVectors
        calc_controlPoints(data[name]) # calculate control points
        calc_panelPoints(data[name])
    
        

    
    return data

############################################################################    

class dataClass(object):
    '''Generates an empty class    
    '''
    def __init__(self):
        pass

class Geometry(object):
    '''Geometry Class. Defines the name, x and y coordinates
    
    Geometry(object) -> object.name, object.x, object.y
    
    '''
    def __init__(self,name,x,y):
        self.name = name
        self.x = x
        self.y = y

def calc_unitVectors(object):
    '''Calculates the unit vectors of the geometry (object)\n
    
    calc_unitVectors(object) -> object.unitNorm, object.unitTang
    '''
    
    norm = normVec((object.points.x[:-1],object.points.y[:-1]),(object.points.x[1:],object.points.y[1:]))
    tang = tangVec((object.points.x[:-1],object.points.y[:-1]),(object.points.x[1:],object.points.y[1:]))
    
    object.unitVectors.norm = Geometry('normal vectors',norm[0],norm[1]) 
    object.unitVectors.tang = Geometry('normal vectors',tang[0],tang[1])
    
    return object
    
def calc_controlPoints(object):
    '''Generates the control points\n
    
    calc_controlPoints(object) -> object.controlPoints.x, object.controlPoints.y
    '''
   
    cpx = (object.points.x[1:] + object.points.x[:-1])/2 + 100*object.unitVectors.norm.x*finfo(float).eps
    cpy = (object.points.y[1:] + object.points.y[:-1])/2 + 100*object.unitVectors.norm.y*finfo(float).eps
   
    object.controlPoints = Geometry('control points',cpx,cpy)
    
    return object
    
def calc_panelPoints(object):
    '''Generates the panel starting and end points
    
    calc_panelPoints(object) -> object.panelStart(.x,.y), object.panelEnd(.x,.y)        
    '''
    
    object.panel.start = Geometry('panel start', object.points.x[:-1], object.points.y[:-1])
    object.panel.end = Geometry('panel end', object.points.x[1:], object.points.y[1:])
    
    return object

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
