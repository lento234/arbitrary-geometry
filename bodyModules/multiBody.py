# -*- coding: utf-8 -*-
"""
Name:       Bodies. [Body and MultiBody modules]
Description:Contains all the data of body and its function
Author:     Lento Manickathan - 1544101
"""

#==============================================================================
# Importing Modules
#==============================================================================

# Standarad scientific module
import numpy as np

# System module: to import
import sys

# Appending path to panelMethod and importing panelMethod
sys.path.append('../panelMethod')

from panelMethod import source2D
from panelMethod import vortex2D

# Custom body modules
from bodyTransform import body2global
import matplotlib.pylab as plt # plotting module
from collocationPoint import collocationPoint # Calculate the location of collocation point
from unitVectors import normalVector, tangentVector # calculate the normal and tangent vector of panels
import time

#==============================================================================
# Multi-body module to store all a multi-body data
#==============================================================================
class multiBody:
    '''
    Multi-body
    '''
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
        
        self.collocationPoint = (self.panelStart + self.panelEnd)/2
        self.vortex_collocationPoint = collocationPoint(self.panelStart,self.panelEnd, -self.normal) # inside the body
        self.source_collocationPoint = collocationPoint(self.panelStart,self.panelEnd, self.normal) # inside the body

    # Plotting function
    def plot(self,attribute,marker='k'):
        
        start = 0
        for num in range(len(self.name)):
            end = start + self.length_geometry[num]

            plt.plot(self.__dict__[attribute][0,start:end],
                     self.__dict__[attribute][1,start:end],marker)
            
            start = end

    def vortexPanel_solve(self, Vinduced, evaluationPoint = 'self'):
        '''
        Solve potential flow using vortex panels
        '''
        
        # Defining collocationPoints: slightly inside
        self.vortex_collocationPoint = collocationPoint(self.panelStart,self.panelEnd, -self.normal) # inside the body
        
        # Calculate RHS - zero tangential
        self.vortex_RHS = vortex2D.RightHandSide(Vinduced,
                                                 self.tangent)
        #start = time.time()
                                         
        # Calculate the influence matrix
        self.vortex_A   = vortex2D.solve(self.vortex_collocationPoint,
                                         self.panelStart,
                                         self.panelEnd,
                                         self.tangent)
                                         
        #print 'vortex:' + str(time.time() - start)
        #start = time.time()
        
        # Solve the vortex Method. Ax=RHS
        self.vortex_gamma = np.linalg.solve(self.vortex_A, self.vortex_RHS)
        
        #print 'vortex:' + str(time.time() - start)
        #start = time.time()

        # To show no transpiration    
        if evaluationPoint is 'self':
            evaluationPoint = self.vortex_collocationPoint
            
        # Calculate induced velocity on body
        self.vortex_V = vortex2D.evaluate(self.vortex_gamma,
                                          evaluationPoint,
                                          self.panelStart,
                                          self.panelEnd)
                                          
       
         
                                          
    #    def vortexNormPanel_solve(self, Vinduced, evaluationPoint = 'self'):
    #        '''
    #        Solve potential flow using vortex panels
    #        '''
    #        
    #        # Defining collocationPoints: slightly inside
    #        self.vortexNorm_collocationPoint = collocationPoint(self.panelStart,self.panelEnd, -self.normal) # inside the body
    #        
    #        # Calculate RHS - zero tangential
    #        self.vortexNorm_RHS = vortex2D.RightHandSide(Vinduced,
    #                                                     self.tangent)
    #                                                 
    #        # Calculate the influence matrix
    #        self.vortexNorm_A   = vortex2D.solve(self.vortexNorm_collocationPoint,
    #                                             self.panelStart,
    #                                             self.panelEnd,
    #                                             self.tangent)
    #                                         
    #        # Solve the vortex Method. Ax=RHS
    #        self.vortexNorm_gamma = np.linalg.solve(self.vortexNorm_A, self.vortexNorm_RHS)
    #        
    #        # To show no transpiration    
    #        if evaluationPoint is 'self':
    #            evaluationPoint = self.vortexNorm_collocationPoint
    #            
    #        # Calculate induced velocity on body
    #        self.vortexNorm_V = vortex2D.evaluate(self.vortexNorm_gamma,
    #                                              evaluationPoint,
    #                                              self.panelStart,
    #                                              self.panelEnd)                                          
                                              
    def sourcePanel_solve(self, Vinduced, evaluationPoint='self'):
        '''
        Solve the panel problem
        '''
        
        # Defining collocationPoint
        self.source_collocationPoint = collocationPoint(self.panelStart, self.panelEnd, self.normal) # + normal
        
        # Calculate RHS - Neumann B.C
        self.source_RHS = source2D.RightHandSide(Vinduced, 
                                                 self.normal)
        
        #start = time.time()
        # Calculate the influence matrix
        self.source_A   = source2D.solve(self.source_collocationPoint, 
                                         self.panelStart,
                                         self.panelEnd,
                                         self.normal)
        

        # Solve the panel Method. Equation: Ax = RHS. Solve for x
        self.source_sigma = np.linalg.solve(self.source_A, self.source_RHS)
        
        #print 'sourcePanel:' + str(time.time() - start)

        # To show no transpiration    
        if evaluationPoint is 'self':
            evaluationPoint = self.source_collocationPoint
         
        start = time.time() 
        # Calculate induced velocity on body
        self.source_V = source2D.evaluate(self.source_sigma,
                                          evaluationPoint,
                                          self.panelStart,
                                          self.panelEnd)
        print str(time.time() - start)

                                          
    def sourceVortexPanel_solve(self, Vinduced, evaluationPoint='self'):
        '''
        '''

        # Defining collocation Point        
        self.source_collocationPoint = collocationPoint(self.panelStart, self.panelEnd, self.normal)
        self.vortex_collocationPoint = collocationPoint(self.panelStart, self.panelEnd, -self.normal)
        
        # Number of collocation points
        n = np.shape(self.collocationPoint)[1]
        
        # Calculate the RHS - Neumann B.C [normal and tangent]
        self.sourceVortex_RHS_normal   = source2D.RightHandSide(Vinduced, self.normal)
        self.sourceVortex_RHS_tangent  = vortex2D.RightHandSide(Vinduced, self.tangent)
        #   Total RHS
        self.sourceVortex_RHS          = np.concatenate([self.sourceVortex_RHS_normal, self.sourceVortex_RHS_tangent])
         
        # Calculate the influence Matrix
        #   normal
        self.sourceVortex_A_normal = source2D.solve(self.source_collocationPoint, self.panelStart, self.panelEnd, self.normal)
        self.sourceVortex_B_normal = vortex2D.solve(self.vortex_collocationPoint, self.panelStart, self.panelEnd, self.normal)
        #   tangent
        self.sourceVortex_A_tangent = source2D.solve(self.source_collocationPoint, self.panelStart, self.panelEnd, self.tangent)
        self.sourceVortex_B_tangent = vortex2D.solve(self.vortex_collocationPoint, self.panelStart, self.panelEnd, self.tangent)
        
        # Combined influence matrix
        self.sourceVortex_C = np.zeros((n*2,n*2)) #[[1,2],[3,4]]
        
        self.sourceVortex_C[:n,:n] = self.sourceVortex_A_normal # 1
        self.sourceVortex_C[:n,n:] = self.sourceVortex_B_normal # 2
        self.sourceVortex_C[n:,:n] = self.sourceVortex_A_tangent # 3 
        self.sourceVortex_C[n:,n:] = self.sourceVortex_B_tangent # 4
    
        #start = time.time()

        # Solve the panel Menthod, combined B.Cs
        self.sourceVortex_SigmaGamma = np.linalg.solve(self.sourceVortex_C, self.sourceVortex_RHS)
        
        #print 'sourceVortexPanel:' + str(time.time() - start)

        #self.sourceVortex_SigmaGamma = np.linalg.lstsq(self.sourceVortex_C, self.sourceVortex_RHS)
        
        self.sourceVortex_sigma = self.sourceVortex_SigmaGamma[:n]
        self.sourceVortex_gamma = self.sourceVortex_SigmaGamma[n:]
        
        # To show no transpiration    
        if evaluationPoint is 'self':
            evaluationPoint_source = self.source_collocationPoint
            evaluationPoint_vortex = self.vortex_collocationPoint
            
        else:
            evaluationPoint_source = evaluationPoint
            evaluationPoint_vortex = evaluationPoint
            
        # Calculate induced velocity on body
        self.sourceVortex_Vsorc = source2D.evaluate(self.sourceVortex_sigma, evaluationPoint_source, self.panelStart, self.panelEnd)        
        self.sourceVortex_Vvort = vortex2D.evaluate(self.sourceVortex_gamma, evaluationPoint_vortex, self.panelStart, self.panelEnd)
        self.sourceVortex_V = self.sourceVortex_Vsorc + self.sourceVortex_Vvort
        

                      
"""                                
    def vortexPanel_solve(self, evaluationPoint = 'self', Velocity_Induced):
        '''
        Solve potential flow using vortex panels
        '''
        
        # Defining collocationPoints: slightly inside
        #self.collocationPoint = collocationPoint(self.panelStart,self.panelEnd, -self.normal) # inside the body
        
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
        self.vortex_gamma = np.linalg.solve(self.vortex_A, self.vortex_RHS)
        
        # To show no transpiration    
        if evaluationPoints is 'self':
            evaluationPoints = self.collocationPoint
            
        # Calculate induced velocity on body
        self.Vinduced = vortex2D.evaluate(self.vortex_gamma,
                                          evaluationPoints,
                                          self.panelStart,
                                          self.panelEnd)
        
"""  
#    def vortexPanel_solve(self, evaluationPoints = 'self', vortexPoints=None, freestream=[[0.],[0.]]):
#        '''
#        Solve potential flow using vortex panels
#        '''
#        
#        # Defining collocationPoints: slightly inside
#        self.collocationPoint = collocationPoint(self.panelStart,self.panelEnd, -self.normal) # inside the body
#        
#        # Calculate RHS - zero tangential
#        self.vortex_RHS = vortex2D.RightHandSide(self.collocationPoint,
#                                                 self.tangent,
#                                                 vortexPoints,
#                                                 freestream)
#                                                 
#        # Calculate the influence matrix
#        self.vortex_A   = vortex2D.solve(self.collocationPoint,
#                                         self.panelStart,
#                                         self.panelEnd,
#                                         self.tangent)
#                                         
#        # Solve the vortex Method. Ax=RHS
#        self.vortex_gamma = np.linalg.solve(self.vortex_A, self.vortex_RHS)
#        
#        # To show no transpiration    
#        if evaluationPoints is 'self':
#            evaluationPoints = self.collocationPoint
#            
#        # Calculate induced velocity on body
#        self.Vinduced = vortex2D.evaluate(self.vortex_gamma,
#                                          evaluationPoints,
#                                          self.panelStart,
#                                          self.panelEnd)
