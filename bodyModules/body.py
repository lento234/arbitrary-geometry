# -*- coding: utf-8 -*-
"""
Name:           body
Description:    Contains all the data of body and its function
Author:         Lento Manickathan - 1544101 - lento.manickathan@gmail.com
"""

#==============================================================================
# Importing modules
#==============================================================================

# Standard scientific module
import numpy as np
import matplotlib.pylab as plt # plotting module

# Custom modules
#from collocationPoint import collocationPoint # Calculate the location of collocation point
#from unitVectors import normalVector, tangentVector # calculate the normal and tangent vector of panels
from bodyTransform import body2global # transform the body coordinates to global geometry
from bodyModules.getCoordinates import getCoordinates # Extracts coordinates

#from bodyModules.generate2DCylinder import generate2DCylinder # Create Circle
import generateCoordinates

#==============================================================================
# Body class: Stores all the data and function related to a given body
#==============================================================================
class body:
    '''
    Stores the input data of the body and calculates more data.
    
    Input
    -----
    name         
        - names of the body in order of input.
    
    shape           
        - tuple: (bodyType, bodyAttribute).
        
            1) bodyType = 'cylinder', bodyAttribute = <number of panels>
                ex: ('cylinder', 100)
                if bodyType is a cylinder, a normalized cylinder coordinates 
                is generated for the given number of panels.
                
            2) bodyType = '<dir>/<airfoilName>.<extension>', bodyAttribute = <if points are clockwise>
                ex: ('geometries/NACA0012.txt',False)
                if bodyType is a dir and filename, then coordinates are 
                imported from the given file path. If bobyAttribute is false, 
                then coordinates are flipped to get clockwise definition of 
                the coordinates.
                              
    chord          
        - chord length or the diameter of the body. Length in meters. [m]
                      
    local_pitch     
        - local pitch (theta) of the airfoil [deg]. Pitch angle should be in 
        degrees and the body is rotated around the 'pivot point'.
                      
    pivot_point     
        - x,y coordinate of center of rotation of the airfoil. The pivot point 
        is considered to be the center of the body. Non-dimensionalized 
        parameter, with respect to chord length. Leading edge (x=0,y=0) and 
        trailing edge: (x=1,y=0).
    
    Returns
    -------
    'object' containing the all the body parameters.
    
    name
        - names of the body in order of input.

    chord
        - chord length of the body. Prefered in [m].
    
    local_pitch     
        - local pitch of the body. In [degrees]. The rotation of the body is 
        about the pivot point of the body.
                      
    pivot_point
        - x,y coordinate of the point of local rotation. The pivot point is 
        considered to be the center of the body. 
    
    geometry_normalized 
        - normalized x,y coordinates of the body in 2D array, row 0 and row 1 
        respectively.
                              
    geometry
        -  x,y-coordinates of the body dimensionalized with chord length, 
        body centered to pivot point and pitch using the local pitch angle 
        around the pivot point. x,y coordinates in row 0 and row 1 respectively.
                                                          
    normal
        - normal x,y component of the body in a 2D array. Calculated using the 
        'normalVector' module from 'geometryModules'.
    
    tangent
         - tangential x-component and y-component of the body. Calculated 
         using the 'tangentVector' module from 'geometryModules'
    
    collocationPoint 
        - collocation points of the geometry in 2D array. Points used to 
        evaluate the body induction terms (source term or vortex term). 
        Collocation points are at the mid of the body panel points. 
        Calculated using 'collocationPoint' module from 'geometryModules'.
    '''
    
    def __init__(self, name, shape, chord=1., local_pitch=0., pivot_point=[0., 0.]):
        # Storing body Parameters
        self.name        = name    # name
        self.chord       = chord  # chord length
        self.local_pitch = local_pitch # pitch angle
        self.pivot_point = np.array([pivot_point]).transpose()  # origin of body w.r.t global geometry
        
        # Generating coordintes
        if shape[0] is 'cylinder':
            self.geometry_normalized = generateCoordinates.cylinder2D(shape[1])
        elif shape[0] is 'block':
            # generate block coordinates
            self.geometry_normalized = generateCoordinates.block2D(shape[1])
        else:
            # body is airfoil; # extract data from file containing normalized coordinates
            # coordinates are given anti-clockwise
            self.geometry_normalized = getCoordinates(filename  = shape[0],
                                                      clockwise = shape[1])          
        
        # Storing body geometry and its parameters
        #   non-normalized geometry (with pitch and correct center)
        self.geometry = body2global((self.geometry_normalized - self.pivot_point)*self.chord, self.local_pitch)        
        
        # Number of coordinate point
        self.noPoints = np.shape(self.geometry)[1]
        
        # Body points, split to panel points - POINTERS
        #self.panelStart = self.geometry[:,:-1] # starting points
        #self.panelEnd   = self.geometry[:,1:]  # ending points
        
        # Unit Vectors of the body
        #self.normal  = normalVector(self.panelStart,  self.panelEnd)         
        #self.tangent = tangentVector(self.panelStart, self.panelEnd)
            
    # Plotting function
    def plot(self,attribute,marker='k'):
        '''
        Plots the attributes of the body.
        
        Input
        -----
        attribute
            - body parameter that needs to plotted. Ex. 'geometry'
        
        marker
            - marker used for plotting. Default is black (k).
        
        Returns
        -------
        figure with attribute plotted. Figure not configured
        '''
        
        plt.plot(self.__dict__[attribute][0],
                 self.__dict__[attribute][1],
                 marker)
    
    # Quiver plotting function        
    def quiver(self,coordinate,attribute):
        '''
        Quiver plot the attribute of the body.chord
        
        Input
        -----
        coordinates - the coordinates of the quiver plot. Standard coordinates
                      of the body: 'geometry'.
                      
        attribute   - body parameter that needs to plotted. Ex. 'geometry'
        
        Returns
        -------
        quiver plot with attributes plotted at given coordinates. Figure not
        configured.        
        '''
                
        plt.quiver(self.__dict__[coordinate][0], self.__dict__[coordinate][1],
                   self.__dict__[attribute][0],  self.__dict__[attribute][1])
 