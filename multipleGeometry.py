# -*- coding: utf-8 -*-
"""
Case:Arbitrary Geometry
Methodology: Solving any geometry using Panel Method
Description: Contains the multipleGeometry.py modules 
@author: lento
"""

# external modules
from numpy import *

############################################################################

# module to calculate panel problem
def panelMethod(geometries, Uinf=0, Winf=0, data=None):
    ''' Calculates the induced velocites for a given geometries or data
    
    sourceTerm(geometries,Uinf,Winf,data) -> data
    
    Parameters
    ---------
    geometries: dictionary of geometries
        Key:
            name of geometry\n
        Content:
            array of x, y coordinates in row 0 and 1 respectively\n
        
    Uinf, Winf (optional): Freestream velocity, optional else zero.
        Uinf:
            Freestream, in x-axis\n
        Winf:
            Freestream, in y-axis\n
        
    data [optional]: data from .sourceTerm can be provided if already calculated
        data.points
        ...
        ...
        data.Sigma
        
    Returns
    -------
    data: 'data' class containing all the data (row 0 = x-axis, row 1 = y-axis)
        data.points: 
            all the points\n
        data.unitVectors (.norm, .tang):
            the unit vectors of panel\n
        data.controlPoints:
            all the control points (mid-point of panel)\n
        data.panel (.start, .end):
            split the panel points to start and end\n
        data.Sigma:
            source term of the problem\n
        
        &&
        
        Induced Velocity at control point
        
        data.Q:
            induced velocity (row 0 = x-axis, row 1 = y-axis)
        data.Qtang:
            induced velocity, tangential to panel
        data.Qnorm:
            induced velocity, normal to panel
        data.Qres:
            induced velocity, resultant
    '''

    # Determining if sigma needs to be calculated    
    if data is None:
        # no data available, so sourceTerm module used to calculate the terms
        data, controlPoints, panelEnd, panelStart, Sigma, unitNorm, unitTang = sourceTerm(geometries, Uinf, Winf, calc_inducedVelocity=True)
    else:
        # already have data
        controlPoints, panelEnd, panelStart, unitNorm, unitTang = reshapeData(data) # reshaping data
        # reshaping Sigma
        Sigma = array([])
        for num in range(len(data)):
            Sigma = concatenate((Sigma,data[data.keys()[num]].Sigma),axis=0)
    
    # Reshaping the Sigma data
    Sigma = tile(Sigma,[shape(controlPoints.x)[0],1])

    #Induced Velocity at control point
    u,w = sor2D(Sigma, (controlPoints.x,controlPoints.y), (panelStart.x, panelStart.y), (panelEnd.x, panelEnd.y))
    
    Q   = concatenate(([sum(u,axis=1)],[sum(w,axis=1)])) + array([[Uinf],[Winf]]) 
    Qt  = Q[0,:]*unitTang.x[:,0] + Q[1,:]*unitTang.y[:,0] # tangential velocity
    Qn  = Q[0,:]*unitNorm.x[:,0] + Q[1,:]*unitNorm.y[:,0] # normal velocity
    Qres = sqrt(pow(Q[0,:],2) + pow(Q[1,:],2)) # resultant velocity
    
    # Storing the data to data[].Q, .Qtang, .Qnorm, .Qres
    start = 0
    for num in range(len(data)):
        end = shape(data[data.keys()[num]].controlPoints)[1]
        data[data.keys()[num]].__setattr__('Q', Q[:,start:(start+end)])
        data[data.keys()[num]].__setattr__('Qtang', Qt[start:(start+end)])
        data[data.keys()[num]].__setattr__('Qnorm', Qn[start:(start+end)])
        data[data.keys()[num]].__setattr__('Qres', Qres[start:(start+end)])
        start = end
    
    return data

# module to calculate the source term (only)
def sourceTerm(geometries, Uinf=0, Winf=0):# calc_inducedVelocity = False):
    ''' Calculates the source term for a given geometries
    
    sourceTerm(geometries,Uinf,Winf) -> data
    
    Parameters
    ---------
    geometries: dictionary of geometries
        Key:
            name of geometry\n
        Content:
            array of x, y coordinates in row 0 and 1 respectively\n
        
    Uinf, Winf (optional): Freestream velocity, optional else zero.
        Uinf:
            Freestream, in x-axis\n
        Winf:
            Freestream, in y-axis\n
        
    calc_inducedVelocity [ignore]:
        
    Returns
    -------
    data: 'data' class containing all the data (row 0 = x-axis, row 1 = y-axis)
        data.points: 
            all the points\n
        data.unitVectors (.norm, .tang):
            the unit vectors of panel\n
        data.controlPoints:
            all the control points (mid-point of panel)\n
        data.panel (.start, .end):
            split the panel points to start and end\n
        data.Sigma:
            source term of the problem\n
    '''
    
    # Initially arranging and organizing the data
    data = organizeData(geometries)
    
    # Reshapes the data accordingally
    controlPoints, panelEnd, panelStart, sigma, unitNorm, unitTang = reshapeData2Solve(data)
    
    # Solving the problem
    u,w = sor2D(sigma, (controlPoints.x, controlPoints.y), (panelStart.x, panelStart.y), (panelEnd.x, panelEnd.y))
    A   = u*unitNorm.x + w*unitNorm.y
    RHS = -(Uinf*unitNorm.x[:,0] + Winf*unitNorm.y[:,0])
    
    Sigma = linalg.solve(A,RHS)
    
    # Storing the data to data[].Sigma
    start = 0
    data.__setattr__('Sigma', {})
    for name in data.geometries:
        end = shape(data.controlPoints[name])[1]
        data.Sigma[name] = Sigma[start: (start + end)]
        start = end
        
    return data

# to be written: work in progress
def inducedVelocities(geometries, collocationPointsX='self', collocationPointsY='self', Uinf=0., Winf=0., data=None):
    '''Calculated induced velocities of any points

    '''
    
    # Solve for sigma and other terms to 
    if data == None:
        data = sourceTerm(geometries, Uinf, Winf)
    
    # reshape Data for calculation
    controlPoints, panelEnd, panelStart, Sigma, unitNorm, unitTang = reshapeData2Calc(data, collocationPointsX, collocationPointsY)
  
    # Solving for induction
    u,w = sor2D(Sigma, (controlPoints.x, controlPoints.y), (panelStart.x, panelStart.y), (panelEnd.x,panelEnd.y))
    
    # Summing the induction velocities
    Qx = sum(u,axis=0) + Uinf
    Qy = sum(w,axis=0) + Winf
    Qres = sqrt(pow(Qx,2) + pow(Qy,2))
    
    # Storing the data to data[].Q, .Qtang, .Qnorm, .Qres
    if collocationPointsX == 'self':
        # If collocation point on body itself
        data.__setattr__('Qx', {})
        data.__setattr__('Qy', {})
        data.__setattr__('Qtang', {}) # unique to self solving
        data.__setattr__('Qnorm', {}) # unique to self solving
        data.__setattr__('Qres', {})
        
        Qtang  = Qx*unitTang.x + Qy*unitTang.y # tangential velocity
        Qnorm  = Qx*unitNorm.x + Qy*unitNorm.y # normal velocity
        
        start = 0
        for name in data.geometries:
            end = shape(data.controlPoints[name])[1]
            
            data.Qx[name]       = Qx[start: (start + end)]
            data.Qy[name]       = Qy[start: (start + end)]
            data.Qres[name]     = Qres[start: (start + end)]
            data.Qtang[name]    = Qtang[start: (start + end)]
            data.Qnorm[name]    = Qnorm[start: (start + end)]

            start = end
    
    else:
        data.__setattr__('Qx', Qx)
        data.__setattr__('Qy', Qy)
        data.__setattr__('Qres', Qres)
    
    return data 

############################################################################
#LIST OF DATA MANIPULATION MODULES

# Empty class generation module 
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

def organizeData(object):
    ''' organizeData(object) -> data (.points, .unitVectors, .controlPoints, .panel (.start, .end))
    
    Method to organize the inital data
    '''
    
    # Creating and empty class    
    data = dataClass()
    
    data.__setattr__('length', len(object))
    data.__setattr__('geometries', object.keys())
    
    data.__setattr__('points', object) # Storing the geometry points
    
    data.__setattr__('unitVectors', dataClass()) # to store norm and tang
    data.unitVectors.__setattr__('norm', {})
    data.unitVectors.__setattr__('tang', {})
    
    data.__setattr__('controlPoints', {})
    
    data.__setattr__('panel', dataClass())
    data.panel.__setattr__('start', {})
    data.panel.__setattr__('end', {})
    
    for name in data.geometries:            
        data.panel.start[name], data.panel.end[name] = calc_panelPoints(data.points[name]) # calculate panel points

        data.unitVectors.norm[name] = calc_normVec((data.panel.start[name][0], data.panel.start[name][1]), (data.panel.end[name][0], data.panel.end[name][1]))
        data.unitVectors.tang[name] = calc_tangVec((data.panel.start[name][0], data.panel.start[name][1]), (data.panel.end[name][0], data.panel.end[name][1]))
        data.controlPoints[name] = calc_controlPoints(data.panel.start[name], data.panel.end[name], data.unitVectors.norm[name] ) # calculate control points
        
    return data
    
def reshapeData2Solve(object):
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
    for name in object.geometries:
        controlPoints.x = concatenate((controlPoints.x, object.controlPoints[name][0]),axis=0)
        controlPoints.y = concatenate((controlPoints.y, object.controlPoints[name][1]),axis=0)
    
        panelStart.x = concatenate((panelStart.x, object.panel.start[name][0]),axis=0)    
        panelStart.y = concatenate((panelStart.y, object.panel.start[name][1]),axis=0)

        panelEnd.x = concatenate((panelEnd.x, object.panel.end[name][0]),axis=0)    
        panelEnd.y = concatenate((panelEnd.y, object.panel.end[name][1]),axis=0)
    
        unitNorm.x = concatenate((unitNorm.x, object.unitVectors.norm[name][0]),axis=0)
        unitNorm.y = concatenate((unitNorm.y, object.unitVectors.norm[name][1]),axis=0)
    
        unitTang.x = concatenate((unitTang.x, object.unitVectors.tang[name][0]),axis=0)
        unitTang.y = concatenate((unitTang.y, object.unitVectors.tang[name][1]),axis=0)
        
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
    
    #Calculating the unit sigma
    sigma = ones(shape(controlPoints.x)) # unit source

    return controlPoints, panelEnd, panelStart, sigma, unitNorm, unitTang
    
def reshapeData2Calc(object, collocationPointsX, collocationPointsY):
    ''' reshapeData2Calc(object, collocationPointsX, collocationPointsY) -> 
    
    Reshape the data to calculate the induced velocity    
    '''
    
    controlPoints   = Geometry(array([]),array([]))
    panelStart      = Geometry(array([]),array([]))
    panelEnd        = Geometry(array([]),array([]))
    unitNorm        = Geometry(array([]),array([]))
    unitTang        = Geometry(array([]),array([]))
    Sigma           = array([])
    
    # concatenating all the data to one matrix
    for name in object.geometries:       
        if collocationPointsX == 'self':            
            controlPoints.x = concatenate((controlPoints.x, object.controlPoints[name][0]), axis=0)
            controlPoints.y = concatenate((controlPoints.y, object.controlPoints[name][1]), axis=0)
            
        unitNorm.x = concatenate((unitNorm.x, object.unitVectors.norm[name][0]), axis=0)
        unitNorm.y = concatenate((unitNorm.y, object.unitVectors.norm[name][1]), axis=0)
    
        unitTang.x = concatenate((unitTang.x, object.unitVectors.tang[name][0]),axis=0)
        unitTang.y = concatenate((unitTang.y, object.unitVectors.tang[name][1]),axis=0)
    
        panelStart.x = concatenate((panelStart.x, object.panel.start[name][0]), axis=0)    
        panelStart.y = concatenate((panelStart.y, object.panel.start[name][1]), axis=0)

        panelEnd.x = concatenate((panelEnd.x, object.panel.end[name][0]), axis=0)    
        panelEnd.y = concatenate((panelEnd.y, object.panel.end[name][1]), axis=0)
    
        Sigma = concatenate((Sigma, object.Sigma[name]), axis=0)
        
    #N = len(controlPoints.x) # total number of control points
    M = len(panelStart.x) # total number of panel points
        
    # Generating points
    if collocationPointsX == 'self':
        controlPoints.x = tile(controlPoints.x, [M,1,1])
        controlPoints.y = tile(controlPoints.y, [M,1,1])
        
        unitNorm.x = tile(unitNorm.x, [M,1,1])
        unitNorm.y = tile(unitNorm.y, [M,1,1])
    
        unitTang.x = tile(unitTang.x, [M,1,1])
        unitTang.y = tile(unitTang.y, [M,1,1])
        
    else:
        controlPoints.x = tile(collocationPointsX, [M,1,1])
        controlPoints.y = tile(collocationPointsY, [M,1,1])
    
    cpShape = shape(controlPoints.x) # Shape of the collocation point

    panelStart.x = tile(reshape(panelStart.x, (M,1,1)), [1,cpShape[1], cpShape[2]])
    panelStart.y = tile(reshape(panelStart.y, (M,1,1)), [1,cpShape[1], cpShape[2]])
    
    panelEnd.x = tile(reshape(panelEnd.x, (M,1,1)), [1,cpShape[1], cpShape[2]])
    panelEnd.y = tile(reshape(panelEnd.y, (M,1,1)), [1,cpShape[1], cpShape[2]])
        
    Sigma = tile(reshape(Sigma, (M,1,1)), [1, cpShape[1], cpShape[2]])
    
    return controlPoints, panelEnd, panelStart, Sigma, unitNorm, unitTang
    
    
############################################################################
# LIST OF CALCULATION MODULES

# module to calculate the source term
def sor2D(sigma, controlPoint, panelStart, panelEnd):
    ''' sor2D(sigma, controlPoints, panelStart, panelEnd) -> u, w
    
    Description: Constant Strength Source Method
    
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

def calc_panelPoints(points):
    ''' calc_panelPoints(points) -> panel start, panel end
    
    Generates the panel starting and end points
    '''
    
    return points[:,:-1], points[:,1:]
    
def calc_controlPoints(start, end, norm):
    ''' calc_controlPoints(points, norm) -> control points
    
    Generates the control points
    '''

    # control points pushed out by delta eps
    cp = add((start + end)/2, 100*norm*finfo(float).eps) 

    return cp
       
# calculate the normal vector
def calc_normVec(panelStart, panelEnd):
    ''' calc_normVec(start, end) -> norm
    
    Description: Finds the normal vector of a line
    
    Parameters: starting point,end point 
      array containing x,y coordinates in row 1 and 2 respectively
      
      Shape: 2 by n
      
    Returns: normal vector of line
      array containing normal unit vector,
      x and y componenent in row 1 and 2 respectively
      
      Shape: 2 by n
      
    '''
    
    # Introducing parameters
    x1, y1 = panelStart
    x2, y2 = panelEnd
    
    r = sqrt((x2-x1)**2 + (y2-y1)**2) #hypotenuse
    
    sinAlpha = (y2-y1)/r
    cosAlpha = (x2-x1)/r
    
    norm = concatenate(([-sinAlpha], [cosAlpha]), axis=0)

    return norm

# calculate the tangential vector
def calc_tangVec(panelStart, panelEnd):
    ''' calc_tangVec(panelStart, panelEnd) -> tang
    
    Description: Finds the tangential vector of a line
    
    Parameters: starting point,end point
      array containing x,y coordinates in row 1 and 2 respectively

      Shape: 2 by n

    Returns: tangent vector of line
      array containing tangent unit vector,
      x and y componenent in row 1 and 2 respectively

      Shape: 2 by n
    '''
    
    # Introducing parameters
    x1, y1 = panelStart
    x2, y2 = panelEnd
    
    r = sqrt((x2 - x1)**2 + (y2 - y1)**2) #hypothenuse

    cosAlpha = (x2 - x1)/r
    sinAlpha = (y2 - y1)/r

    tang = concatenate(([cosAlpha], [sinAlpha]))
    
    return tang
    
############################################################################
# TRANSFORMATION MODULES

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
    
    r = sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    cosAlpha = (x2 - x1)/r
    sinAlpha = (y2 - y1)/r
    
    xp = cosAlpha*(xG - x1)  + sinAlpha*(yG - y1)
    yp = -sinAlpha*(xG - x1) + cosAlpha*(yG - y1)
    
    return xp, yp
    
# to transfer velocity from panel to global geometry
def panel2global(velocityPanel, origin, end):
    '''panel2global(velocityPanel, origin, end) -> u, w
    
    Transforms panel velocity component to global direction 
    
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
    up, wp = velocityPanel # can be coordinate or velocity component(if later x=u, y=w)
    x1, y1 = origin
    x2, y2 = end
    
    #Figure 11.17
    r = sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    cosAlpha = (x2 - x1)/r
    sinAlpha = (y2 - y1)/r
    
    #Equation 10.7
    u = cosAlpha*up - sinAlpha*wp
    w = sinAlpha*up + cosAlpha*wp
    
    return u, w
