# -*- coding: utf-8 -*-
"""
Name:       Bodies. [Body and MultiBody modules]
Description:Contains all the data of body and its function
Author:     Lento Manickathan - 1544101
"""

import numpy as np

from bodyTransform import body2global
from pylab import plot

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
            self.length_geometry.append(np.shape(data['body'].geometry)[1]) # number of geometry points
            self.length_panel.append(np.shape(data['body'].collocationPoint)[1]) # number of panel points
            
            
        # Initializing the arrays for storing all the data [reduce computational load]
        self.chord       = np.zeros(len(args)) # 1-D array of chords
        self.local_pitch = np.zeros(len(args))
        self.global_pitch= np.zeros(len(args))
        self.location    = np.zeros((2,len(args))) 
        
        self.geometry    = np.zeros((2,sum(self.length_geometry))) # x,y coordinates in row 0 and 1
        self.panelStart  = np.zeros((2,sum(self.length_panel)))
        self.panelEnd    = np.zeros((2,sum(self.length_panel)))
        self.normal      = np.zeros((2,sum(self.length_panel)))
        self.tangent     = np.zeros((2,sum(self.length_panel)))
        self.collocationPoint = np.zeros((2,sum(self.length_panel)))
        
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
            self.panelStart[:,start_panel:end_panel]       = body2global(data['body'].panelStart, self.global_pitch[i], self.location[:,i])
            self.panelEnd[:,start_panel:end_panel]         = body2global(data['body'].panelEnd,   self.global_pitch[i], self.location[:,i])
            self.normal[:,start_panel:end_panel]           = data['body'].normal # unit vector is panel specific
            self.tangent[:,start_panel:end_panel]          = data['body'].tangent # unit vector is panel specific
            self.collocationPoint[:,start_panel:end_panel] = body2global(data['body'].collocationPoint, self.global_pitch[i], self.location[:,i])
            
            # Preparing for the next iteration
            start_geometry = end_geometry
            start_panel = end_panel
            i+=1
            
    # Plotting function
    def plot(self,attribute,marker='k'):
        
        start = 0
        for num in range(len(self.name)):
            end = start + self.length_geometry[num]

            plot(self.__dict__[attribute][0,start:end],
                 self.__dict__[attribute][1,start:end],marker)
            
            start = end
    