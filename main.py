# -*- coding: utf-8 -*-
"""
Case: Arbitrary Geometry\n
Methodology: Solving any geometry using Panel Method\n

"""

from pylab import *
ion()

from numpy import *
#from scipy import *

import multipleGeometry

     
# Control Parameters
Uinf = 10.0
Winf = 0.0

# Dictionary of Geometries
Geometries = {'a': array([[-1.0,-0.75,-0.5,-0.75,-1.0], # x-coordinate
                          [0.0,0.25,0.0,-0.25,0.0]]),# y- coordinate
              'b': array([[0.5,0.75,1.0,0.75,0.5],
                          [0.0,0.25,0.0,-0.25,0.0]])}#,
              #'c': array([[1.5,1.75,2.0,1.75,1.5],
              #            [0.0,0.25,0.0,-0.25,0.0]]),
              #'d': array([[2.5,2.75,3.0,2.75,2.5],
              #            [0.0,0.25,0.0,-0.25,0.0]])}

collocationPointsX, collocationPointsY = meshgrid(linspace(-2,2,50),linspace(-0.5,0.5,20))


# Calculating the data (from multipleGeometry module)
#data = multipleGeometry.sourceTerm(Geometries,Uinf=10.)
#data, Qx, Qy, Qres = multipleGeometry.inducedVelocities(Geometries, Uinf=10.)
data = multipleGeometry.inducedVelocities(Geometries, Uinf=10.)

#data = multipleGeometry.panelMethod(Geometries,Uinf=10.)
#data = multipleGeometry.inducedVelocities(Geometries, collocationPointsX, collocationPointsY, Uinf=10.0)

# Plotting Data
#figure(1)
#for name in data.geometries:
#    plot(data.points[name][0],data.points[name][1],'k')
#    quiver(data.controlPoints[name][0], data.controlPoints[name][1], data.Qx[name],data.Qy[name],data.Qres[name])
##contourf(collocationPointsX, collocationPointsY, Qres)
#xlabel('x-coordinate [-]')
#ylabel('y-coordinate [-]')
#title('Multiple Geometries [%d bodies]' %(data.length))
#axis('equal')
#grid()