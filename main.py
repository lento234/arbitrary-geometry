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
airfoilC = Geometry('airfoilC',array([1.5,1.75,2.0,1.75,1.5]),array([0.0,0.25,0.0,-0.25,0.0]))
airfoilD = Geometry('airfoilD',array([2.5,2.75,3.0,2.75,2.5]),array([0.0,0.25,0.0,-0.25,0.0]))

# aribitrary.py module
geometries = (airfoilA,airfoilB,airfoilC,airfoilD)

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

unitNorm.x = tile(array([unitNorm.x]).transpose(),[1,shape(controlPoints.x)[1]])
unitNorm.y = tile(array([unitNorm.y]).transpose(),[1,shape(controlPoints.y)[1]])

unitTang.x = tile(array([unitTang.x]).transpose(),[1,shape(controlPoints.x)[1]])
unitTang.y = tile(array([unitTang.y]).transpose(),[1,shape(controlPoints.x)[1]])

sigma = ones([shape(controlPoints.x)[0],shape(panelStart.x)[1]])

####################################################################
'''             Solving the problem                             '''

u,w = sor2D(sigma, (controlPoints.x,controlPoints.y), (panelStart.x, panelStart.y), (panelEnd.x, panelEnd.y))
A = u*unitNorm.x + w*unitNorm.y
RHS = -(Uinf*unitNorm.x[:,0] + Winf*unitNorm.y[:,0])

Sigma = solve(A,RHS)

####################################################################
'''               Calculating the results                        '''

Sigma = tile(Sigma,[shape(controlPoints.x)[0],1])
#Induced Velocity at control point
u,w = sor2D(Sigma, (controlPoints.x,controlPoints.y), (panelStart.x, panelStart.y), (panelEnd.x, panelEnd.y))
Q     = concatenate(([sum(u,axis=1)],[sum(w,axis=1)])) + array([[Uinf],[Winf]]) 
Qt    = Q[0,:]*unitTang.x[:,0] + Q[1,:]*unitTang.y[:,0] # tangential velocity
Qn    = Q[0,:]*unitNorm.x[:,0] + Q[1,:]*unitNorm.y[:,0] # normal velocity
Qres  = sqrt(pow(Q[0,:],2) + pow(Q[1,:],2)) # resultant velocity

#Storing Data
class finalResult(object):
    def __init__(self,name):
        self.name = name
        
result = finalResult('Results')
result.__setattr__('Geometries',array([add(range(len(geometries)),1)]).transpose())

result.__setattr__('controlPoints',finalResult('Panel points'))
result.controlPoints.__setattr__('x',reshape(controlPoints.x[:,0],(len(geometries),len(geometries[0].cp.x))))
result.controlPoints.__setattr__('y',reshape(controlPoints.y[:,0],(len(geometries),len(geometries[0].cp.y))))

result.__setattr__('Q',finalResult('Induced velocity'))
result.Q.__setattr__('x',reshape(Q[0],(len(geometries),len(geometries[0].cp.x))))
result.Q.__setattr__('y',reshape(Q[1],(len(geometries),len(geometries[0].cp.y))))

result.__setattr__('Qt',reshape(Qt,(len(geometries),len(geometries[0].cp.x))))
result.__setattr__('Qn',reshape(Qn,(len(geometries),len(geometries[0].cp.x))))
result.__setattr__('Qres',reshape(Qres,(len(geometries),len(geometries[0].cp.x))))


figure(1)
plot(result.controlPoints.x,result.controlPoints.y,'b.')
quiver(result.controlPoints.x,result.controlPoints.y,result.Q.x,result.Q.y,result.Qres)
xlabel('x-coordinate [-]')
ylabel('y-coordinate [-]')
axis([-2, 2, -2, 2])
axis('equal')
grid()

