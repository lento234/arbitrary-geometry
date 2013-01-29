# -*- coding: utf-8 -*-
"""
Case: Arbitrary Geometry\n
Methodology: Solving any geometry using Panel Method\n

"""

from pylab import *
ion()

from numpy import *
from scipy import *

# Custom Modules
from Vort2D import *
from unitVec import *

# Geometry Class
class Geometry(object):
    '''
    Geometry Class. Defines the name, x and y coordinates
    
    '''
    def __init__(self,name,x,y):
        self.name = name
        self.x = x
        self.y = y
        
    def unitVectors(self):
        self.unitNorm = normVec((self.x[:-1],self.y[:-1]),(self.x[1:],self.y[1:]))
        self.unitTang = tangVec((self.x[:-1],self.y[:-1]),(self.x[1:],self.y[1:]))
        
    def controlPoints(self):   
        '''
        Generates the control points
        
        '''
        cpx = (self.x[1:]+self.x[:-1])/2 + 100*self.unitNorm[0]*finfo(float).eps
        cpy = (self.y[1:]+self.y[:-1])/2 + 100*self.unitNorm[1]*finfo(float).eps
        self.cp = Geometry('control points',cpx,cpy)
    def panelPoints(self):
        '''
        Generates the panel starting and end points
        
        '''
        self.panelStart = Geometry('panel start',self.x[:-1],self.y[:-1])
        self.panelEnd = Geometry('panel end',self.x[1:],self.y[1:])    
         
        
# Control Parameters
Uinf = 10.0
Winf = 0.0

# Defining Geometry      
airfoilA = Geometry('airfoilA',array([-1.0,-0.75,-0.5,-0.75,-1.0]),array([0.0,0.25,0.0,-0.25,0.0]))
airfoilB = Geometry('airfoilB',array([0.5,0.75,1.0,0.75,0.5]),array([0.0,0.25,0.0,-0.25,0.0]))


# aribitrary.py module
geometries = (airfoilA,airfoilB)

controlPoints = Geometry('control point',array([]),array([]))
panelStart = Geometry('panel start',array([]),array([]))
panelEnd = Geometry('panel end',array([]),array([]))
unitNorm = Geometry('normal vector',array([]),array([]))
unitTang = Geometry('tang vector',array([]),array([]))

for num in range(len(geometries)):
    
    # Calling the methods to calculate control points and panel edge points
    geometries[num].unitVectors()
    geometries[num].controlPoints()
    geometries[num].panelPoints()
    
    # Concatenating the geometries
    controlPoints.x = concatenate((controlPoints.x,geometries[num].cp.x),axis=0)
    controlPoints.y = concatenate((controlPoints.y,geometries[num].cp.y),axis=0)
    
    panelStart.x = concatenate((panelStart.x,geometries[num].panelStart.x),axis=0)    
    panelStart.y = concatenate((panelStart.y,geometries[num].panelStart.y),axis=0)
    
    panelEnd.x = concatenate((panelEnd.x,geometries[num].panelEnd.x),axis=0)    
    panelEnd.y = concatenate((panelEnd.y,geometries[num].panelEnd.y),axis=0)
    
    unitNorm.x = concatenate((unitNorm.x,geometries[num].unitNorm[0]),axis=0)
    unitNorm.y = concatenate((unitNorm.y,geometries[num].unitNorm[1]),axis=0)
    unitTang.x = concatenate((unitTang.x,geometries[num].unitTang[0]),axis=0)
    unitTang.y = concatenate((unitTang.y,geometries[num].unitTang[1]),axis=0)
    
    
# Meshing the data
controlPoints.x = tile(array([controlPoints.x]).transpose(),[1,len(panelStart.x)])
controlPoints.y = tile(array([controlPoints.y]).transpose(),[1,len(panelStart.y)])

panelStart.x = tile(panelStart.x,[shape(controlPoints.x)[0],1])
panelStart.y = tile(panelStart.y,[shape(controlPoints.y)[0],1])

panelEnd.x = tile(panelEnd.x,[shape(controlPoints.x)[0],1])
panelEnd.y = tile(panelEnd.y,[shape(controlPoints.y)[0],1])

unitNorm.x = tile(array([unitNorm.x]).transpose(),[1,len(panelStart.x)])
unitNorm.y = tile(array([unitNorm.y]).transpose(),[1,len(panelStart.y)])

unitTang.x = tile(array([unitTang.x]).transpose(),[1,len(panelStart.x)])
unitTang.y = tile(array([unitTang.y]).transpose(),[1,len(panelStart.y)])

sigma = ones([shape(controlPoints.x)[0],shape(panelStart.x)[1]])

# Solving the problem
u,w = sor2D(sigma, (controlPoints.x,controlPoints.y), (panelStart.x, panelStart.y), (panelEnd.x, panelEnd.y))
A = u*unitNorm.x + w*unitNorm.y
RHS = -(Uinf*unitNorm.x[:,0] + Winf*unitNorm.y[:,0])

Sigma = solve(A,RHS)

#coordinates.x = reshape(coordinates.x,(len(geometries),1,len(geometries[0].x)))
#coordinates.y = reshape(coordinates.y,(len(geometries),1,len(geometries[0].y)))
#coordinates.cpx = reshape(coordinates.cpx,(len(geometries),1,len(geometries[0].cpx)))
#coordinates.cpy = reshape(coordinates.cpy,(len(geometries),1,len(geometries[0].cpy)))



