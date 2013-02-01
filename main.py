# -*- coding: utf-8 -*-
"""
Case: Arbitrary Geometry\n
Methodology: Solving any geometry using Panel Method\n

"""

from pylab import *
ion()

from numpy import *

#from scipy import *

# Custom Modules
#from Vort2D import *
#from unitVec import *

from multipleGeometry import sourceTerm

     
# Control Parameters
Uinf = 10.0
Winf = 0.0

# Dictionary of Geometries
#Geometries = {'a':(array([-1.0,-0.75,-0.5,-0.75,-1.0]),array([0.0,0.25,0.0,-0.25,0.0])),
#              'b':(array([0.5,0.75,1.0,0.75,0.5]),array([0.0,0.25,0.0,-0.25,0.0]))}

Geometries = {'a':array([[-1.0,-0.75,-0.5,-0.75,-1.0], # x-coordinate
                         [0.0,0.25,0.0,-0.25,0.0]]), # y- coordinate
              'b':array([[0.5,0.75,1.0,0.75,0.5],
                         [0.0,0.25,0.0,-0.25,0.0]])}


data = sourceTerm(Geometries)


# Defining Geometry      
#airfoilA = Geometry('airfoilA',array([-1.0,-0.75,-0.5,-0.75,-1.0]),array([0.0,0.25,0.0,-0.25,0.0]))
#airfoilB = Geometry('airfoilB',array([0.5,0.75,1.0,0.75,0.5]),array([0.0,0.25,0.0,-0.25,0.0]))
#airfoilC = Geometry('airfoilC',array([1.5,1.75,2.0,1.75,1.5]),array([0.0,0.25,0.0,-0.25,0.0]))
#airfoilD = Geometry('airfoilD',array([2.5,2.75,3.0,2.75,2.5]),array([0.0,0.25,0.0,-0.25,0.0]))

# aribitrary.py module

#geometries = (airfoilA,airfoilB)#,airfoilC,airfoilD)

#controlPoints = Geometry('control point',array([]),array([]))
#panelStart = Geometry('panel start',array([]),array([]))
#panelEnd = Geometry('panel end',array([]),array([]))
#unitNorm = Geometry('normal vector',array([]),array([]))
#unitTang = Geometry('tang vector',array([]),array([]))

#def calc_unitVectors(object):
#    '''Calculates the unit vectors of the geometry (object)\n
#    
#    calc_unitVectors(object) -> object.unitNorm, object.unitTang
#    '''
#    object.unitNorm = normVec((object.x[:-1],object.y[:-1]),(object.x[1:],object.y[1:]))
#    object.unitTang = tangVec((object.x[:-1],object.y[:-1]),(object.x[1:],object.y[1:]))
#    
#    return object
#    
#def calc_controlPoints(object):
#    '''Generates the control points\n
#    
#    calc_controlPoints(object) -> object.controlPoints.x, object.controlPoints.y
#    '''
#   
#    cpx = (object.x[1:] + object.x[:-1])/2 + 100*object.unitNorm[0]*finfo(float).eps
#    cpy = (object.y[1:] + object.y[:-1])/2 + 100*object.unitNorm[1]*finfo(float).eps
#   
#    object.controlPoints = Geometry('control points',cpx,cpy)
#    
#def calc_panelPoints(object):
#    '''Generates the panel starting and end points
#    
#    
#        
#    '''
#    
#    object.panelStart = Geometry('panel start', object.x[:-1], object.y[:-1])
#    object.panelEnd = Geometry('panel end', object.x[1:], object.y[1:])
    

#for num in range(len(geometries)):
#    
#    # Calling the methods to calculate control points and panel edge points
#    calc_unitVectors(geometries[num])
#    calc_controlPoints(geometries[num])
#    calc_panelPoints(geometries[num])
#    
#    # Concatenating the geometries
#    controlPoints.x = concatenate((controlPoints.x,geometries[num].controlPoints.x),axis=0)
#    controlPoints.y = concatenate((controlPoints.y,geometries[num].controlPoints.y),axis=0)
#    
#    panelStart.x = concatenate((panelStart.x,geometries[num].panelStart.x),axis=0)    
#    panelStart.y = concatenate((panelStart.y,geometries[num].panelStart.y),axis=0)
#    
#    panelEnd.x = concatenate((panelEnd.x,geometries[num].panelEnd.x),axis=0)    
#    panelEnd.y = concatenate((panelEnd.y,geometries[num].panelEnd.y),axis=0)
#    
#    unitNorm.x = concatenate((unitNorm.x,geometries[num].unitNorm[0]),axis=0)
#    unitNorm.y = concatenate((unitNorm.y,geometries[num].unitNorm[1]),axis=0)
#    unitTang.x = concatenate((unitTang.x,geometries[num].unitTang[0]),axis=0)
#    unitTang.y = concatenate((unitTang.y,geometries[num].unitTang[1]),axis=0)
#    
#    
## Meshing the data
#controlPoints.x = tile(array([controlPoints.x]).transpose(),[1,len(panelStart.x)])
#controlPoints.y = tile(array([controlPoints.y]).transpose(),[1,len(panelStart.y)])
#
#panelStart.x = tile(panelStart.x,[shape(controlPoints.x)[0],1])
#panelStart.y = tile(panelStart.y,[shape(controlPoints.y)[0],1])
#
#panelEnd.x = tile(panelEnd.x,[shape(controlPoints.x)[0],1])
#panelEnd.y = tile(panelEnd.y,[shape(controlPoints.y)[0],1])
#
#unitNorm.x = tile(array([unitNorm.x]).transpose(),[1,shape(controlPoints.x)[1]])
#unitNorm.y = tile(array([unitNorm.y]).transpose(),[1,shape(controlPoints.y)[1]])
#
#unitTang.x = tile(array([unitTang.x]).transpose(),[1,shape(controlPoints.x)[1]])
#unitTang.y = tile(array([unitTang.y]).transpose(),[1,shape(controlPoints.x)[1]])
#
#sigma = ones([shape(controlPoints.x)[0],shape(panelStart.x)[1]])
#
#####################################################################
#'''             Solving the problem                             '''
#
#u,w = sor2D(sigma, (controlPoints.x,controlPoints.y), (panelStart.x, panelStart.y), (panelEnd.x, panelEnd.y))
#A = u*unitNorm.x + w*unitNorm.y
#RHS = -(Uinf*unitNorm.x[:,0] + Winf*unitNorm.y[:,0])
#
#Sigma = solve(A,RHS)
#
#####################################################################
#'''               Calculating the results                        '''
#
#Sigma = tile(Sigma,[shape(controlPoints.x)[0],1])
##Induced Velocity at control point
#u,w = sor2D(Sigma, (controlPoints.x,controlPoints.y), (panelStart.x, panelStart.y), (panelEnd.x, panelEnd.y))
#Q     = concatenate(([sum(u,axis=1)],[sum(w,axis=1)])) + array([[Uinf],[Winf]]) 
#Qt    = Q[0,:]*unitTang.x[:,0] + Q[1,:]*unitTang.y[:,0] # tangential velocity
#Qn    = Q[0,:]*unitNorm.x[:,0] + Q[1,:]*unitNorm.y[:,0] # normal velocity
#Qres  = sqrt(pow(Q[0,:],2) + pow(Q[1,:],2)) # resultant velocity
#
##Storing Data
#class finalResult(object):
#    def __init__(self,name):
#        self.name = name
#        
#result = finalResult('Results')
#result.__setattr__('Geometries',array([add(range(len(geometries)),1)]).transpose())
#
#result.__setattr__('controlPoints',finalResult('Panel points'))
#result.controlPoints.__setattr__('x',reshape(controlPoints.x[:,0],(len(geometries),len(geometries[0].controlPoints.x))))
#result.controlPoints.__setattr__('y',reshape(controlPoints.y[:,0],(len(geometries),len(geometries[0].controlPoints.y))))
#
#result.__setattr__('Q',finalResult('Induced velocity'))
#result.Q.__setattr__('x',reshape(Q[0],(len(geometries),len(geometries[0].controlPoints.x))))
#result.Q.__setattr__('y',reshape(Q[1],(len(geometries),len(geometries[0].controlPoints.y))))
#
#result.__setattr__('Qt',reshape(Qt,(len(geometries),len(geometries[0].controlPoints.x))))
#result.__setattr__('Qn',reshape(Qn,(len(geometries),len(geometries[0].controlPoints.x))))
#result.__setattr__('Qres',reshape(Qres,(len(geometries),len(geometries[0].controlPoints.x))))
#
#
#figure(1)
#plot(result.controlPoints.x,result.controlPoints.y,'b.')
#quiver(result.controlPoints.x,result.controlPoints.y,result.Q.x,result.Q.y,result.Qres)
#xlabel('x-coordinate [-]')
#ylabel('y-coordinate [-]')
#axis([-2, 2, -2, 2])
#axis('equal')
#grid()

