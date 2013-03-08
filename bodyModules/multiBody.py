# -*- coding: utf-8 -*-
"""
Name:       Bodies. [Body and MultiBody modules]
Description:Contains all the data of body and its function
Author:     Lento Manickathan - 1544101
"""

import numpy as np
import sys

sys.path.append('../panelMethod')

from bodyTransform import body2global
import matplotlib.pylab as plt # plotting module
from collocationPoint import collocationPoint # Calculate the location of collocation point
from unitVectors import normalVector, tangentVector # calculate the normal and tangent vector of panels

from panelMethod import source2D
from panelMethod import vortex2D
class multiBody:
    def __init__(self,*args):
        
        # Extracting initial data
        # Initializing the lists of variables
        self.name = [] 
        self.length_geometry = []  
        self.length_panel    = []
        # Scanning through the 'args' and getting the data
        for data in args:
            # Extracting the global parameters of the body
            
            # Extracting initial data from the bodies
            self.name.append(data['body'].name) # body names
            self.length_geometry.append(data['body'].noPoints) # number of geometry points
            self.length_panel.append(data['body'].noPoints-1) # number of panel points
            
        # Initializing the arrays for storing all the data [reduce computational load]
        self.chord       = np.zeros(len(args)) # 1-D array of chords
        self.local_pitch = np.zeros(len(args))
        self.global_pitch= np.zeros(len(args))
        self.location    = np.zeros((2,len(args))) 
        
        self.geometry    = np.zeros((2,sum(self.length_geometry))) # x,y coordinates in row 0 and 1
        self.panelStart  = np.zeros((2,sum(self.length_panel)))
        self.panelEnd    = np.zeros((2,sum(self.length_panel)))
        
        # Storing all the data
        
        # Pointer parameters
        i = 0
        start_geometry = 0
        start_panel = 0
        # Looping through the available data
        for data in args:
            
            # Calculating the shape of 'sub array'. Finding the end point
            end_geometry = start_geometry + self.length_geometry[i]
            end_panel = start_panel + self.length_panel[i]
            
            # Storing all the data together into an array
            self.chord[i]       = data['body'].chord
            self.local_pitch[i] = data['body'].local_pitch
            self.global_pitch[i]= data['global_pitch']
            self.location[:,i]  = np.array([data['location']])
            
            # Transforming the points and vectors to global coordinates.          
            self.geometry[:,start_geometry:end_geometry]   = body2global(data['body'].geometry,   self.global_pitch[i], self.location[:,i])
            self.panelStart[:,start_panel:end_panel]       = self.geometry[:,start_geometry:(end_geometry-1)] # starting points
            self.panelEnd[:,start_panel:end_panel]         = self.geometry[:,(1+start_geometry):end_geometry]
              
            # Preparing for the next iteration
            start_geometry = end_geometry
            start_panel = end_panel
            i+=1
        
        self.normal  = normalVector(self.panelStart,  self.panelEnd)
        self.tangent = tangentVector(self.panelStart,  self.panelEnd)
            
    # Plotting function
    def plot(self,attribute,marker='k'):
        
        start = 0
        for num in range(len(self.name)):
            end = start + self.length_geometry[num]

            plt.plot(self.__dict__[attribute][0,start:end],
                     self.__dict__[attribute][1,start:end],marker)
            
            start = end
            
    def sourcePanel_solve(self, evaluationPoints='self', vortexPoints=None, freestream=[0.,0.]):
        '''
        Solve the panel problem
        '''
        
        # Defining collocationPoint
        self.collocationPoint = collocationPoint(self.panelStart, self.panelEnd, self.normal) # + normal
        
        # Calculate RHS - Neumann B.C
        self.source_RHS = source2D.RightHandSide(self.collocationPoint, 
                                                 self.normal, 
                                                 vortexPoints, 
                                                 freestream)
        
        # Calculate the influence matrix
        self.source_A   = source2D.solve(self.collocationPoint, 
                                         self.panelStart,
                                         self.panelEnd,
                                         self.normal)
        
        # Solve the panel Method. Equation: Ax = RHS. Solve for x
        self.source_sigma = np.linalg.solve(self.source_A,self.source_RHS)
    
        # To show no transpiration    
        if evaluationPoints is 'self':
            evaluationPoints = self.collocationPoint
            
        # Calculate induced velocity on body
        self.source_Vinduced = source2D.evaluate(self.source_sigma,
                                                 evaluationPoints,
                                                 self.panelStart,
                                                 self.panelEnd)
        
    def vortexPanel_solve(self, evaluationPoints = 'self', vortexPoints=None, freestream=[[0.],[0.]]):
        '''
        Solve potential flow using vortex panels
        '''
        
        # Defining collocationPoints: slightly inside
        self.collocationPoint = collocationPoint(self.panelStart,self.panelEnd, -self.normal) # inside the body
        
        # Calculate RHS - zero tangential
        self.vortex_RHS = vortex2D.RightHandSide(self.collocationPoint,
                                                 self.tangent,
                                                 vortexPoints,
                                                 freestream)
                                                 
        # Calculate the influence matrix
        self.vortex_A   = vortex2D.solve(self.collocationPoint,
                                         self.panelStart,
                                         self.panelEnd,
                                         self.tangent)
                                         
        # Solve the vortex Method. Ax=RHS
        self.vortex_gamma = np.linalg.solve(self.vortex_A,self.vortex_RHS)
        
        # To show no transpiration    
        if evaluationPoints is 'self':
            evaluationPoints = self.collocationPoint
            
        # Calculate induced velocity on body
        self.Vinduced = vortex2D.evaluate(self.vortex_gamma,
                                          evaluationPoints,
                                          self.panelStart,
                                          self.panelEnd)
        
        
    