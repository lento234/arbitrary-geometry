# -*- coding: utf-8 -*-
"""
Name:       Bodies. [Body and MultiBody modules]
Description:Contains all the data of body and its function
Author:     Lento Manickathan - 1544101
"""

import numpy as np

from collocationPoint import collocationPoint 
from unitVectors import normalVector, tangentVector
from bodyTransform import body2global


class body:
    def __init__(self, name, coordinates, chord=1., local_pitch=0., pivot_point=[0., 0.]):
        '''
        Stores the input data of the body and calculates more data.
        
        Input
        -----
        coordinates     - normalized x,y-coordinates of the body respectively. 
                          1D  array of grid points. [-]
                                                              
        chord           - chord length or the diameter of the body. Length in 
                          meters. [m] 
                          
        theta           - local pitch of the airfoil. pitch angle should be in
                          degrees. [deg] 
                          
        pivot_point     - Center of rotation of the airfoil. If not stated, the
                          center of local pitch will be at the 1/4th of the 
                          chord length. [m]       
        
        Parameters
        ----------
        geometry        - x-coordinate and y-coordinate
        
        chord           - chord length
        
        normal          - normal x-component and y-component of the body. Calculated
                          using the 'normalVector' module from 'geometryModules'.
        
        tangent         - tangential x-component and y-component of the body. Calculated
                          using the 'tangentVector' module from 'geometryModules'
        
        collocationPoint - collocation points of the geometry. Points used to
                          evaluate the body induction terms. Collocation points are 
                          at the mid of the body panel points. Calculated using 
                          'collocationPoint' module from 'geometryModules'.
        '''
        
        # Storing body Parameters
        self.name  = name    # name
        self.chord = chord  # chord length
        self.local_pitch = local_pitch # pitch angle
        self.pivot_point = np.array([pivot_point]).transpose()  # origin of body w.r.t global geometry
        
        self.geometry_normalized = coordinates # normalized coordinates
         
        # Storing body geometry and its parameters
        #self.geometry = body2global(np.multiply(np.add(self.geometry_normalized, -self.pivot_point), self.chord), [0.,0.], self.local_pitch)
        self.geometry = body2global((self.geometry_normalized - self.pivot_point)*self.chord, self.local_pitch)        
        # non-normalized geometry (with pitch and correct center)
        
        # Body points, split to panel points - POINTERS
        self.panelStart = self.geometry[:,:-1] # starting points
        self.panelEnd   = self.geometry[:,1:]  # ending points
        
        # Unit Vectors of the body
        self.normal     = normalVector(self.panelStart,  self.panelEnd)         
        self.tangent    = tangentVector(self.panelStart, self.panelEnd)
        
        # Calculating the collocation points.
        self.collocationPoint = collocationPoint(self.panelStart, self.panelEnd, self.normal)
                
        # print method
        # get method