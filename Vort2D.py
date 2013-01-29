# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 15:45:18 2012
@module: Vort2D
@description: Evaluate the Vorticity in 2D
@input: Vorticity (Gamma), Control point (xCP, yCP),location of induction (x, y)
@output: Induction velocity (u,w)
@author: lento
"""

from numpy import *
#from scipy import *
#from pylab import *

from unitVec import *

# Equation (11.1) of Low speed aerodynamics,  
def sor2D(sigma,controlPoint,panelStart,panelEnd):
    '''Description: Constant Strength Source Method
    
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

def pointVor(Gamma,controlPoint,vortexPoint):
    ''' Two-dimensional point vortex:
        Equation (10.9 & 10.10) '''
    
    # In global coordinates
    x,y     = controlPoint
    x0,y0   = vortexPoint
    
    # Defining preliminary variables
    r       = (x-x0)**2 + (y-y0)**2
    
    # Equation 10.9 and 10.10
    u       = (Gamma/(2*pi))*((y-y0)/r)
    w       = (-Gamma/(2*pi))*((x-x0)/r)
    
    return u,w
    
def pointVorMesh(Gamma,controlPoints,vortexPoint):
    ''' Mesh Evaluation'''
    # In global coordinates
    x,y     = controlPoints
    x0,y0   = vortexPoint
    
    # Defining preliminary variables
    r       = (x-x0)**2 + (y-y0)**2
    # Equation 10.9 and 10.10
    u      = (Gamma/(2*pi))*((y-y0)/r)
    w      = (-Gamma/(2*pi))*((x-x0)/r)
    
    # replacing Nan value
    u[isnan(u)] = 0.
    w[isnan(w)] = 0.
    
    return u,w
