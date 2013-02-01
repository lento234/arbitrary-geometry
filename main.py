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
                          [0.0,0.25,0.0,-0.25,0.0]]), # y- coordinate
              'b': array([[0.5,0.75,1.0,0.75,0.5],
                          [0.0,0.25,0.0,-0.25,0.0]]),
              'c': array([[1.5,1.75,2.0,1.75,1.5],
                          [0.0,0.25,0.0,-0.25,0.0]]),
              'd': array([[2.5,2.75,3.0,2.75,2.5],
                          [0.0,0.25,0.0,-0.25,0.0]])}


# Calculating the data (from multipleGeometry module)
#data = multipleGeometry.sourceTerm(Geometries,Uinf=10.)
data = multipleGeometry.panelMethod(Geometries,Uinf=10.)

# Plotting Data
figure(1)
for num in range(len(data)):
    plot(data[data.keys()[num]].controlPoints[0],data[data.keys()[num]].controlPoints[1],'b.')
    plot(data[data.keys()[num]].points[0],data[data.keys()[num]].points[1],'g')
    quiver(data[data.keys()[num]].controlPoints[0],data[data.keys()[num]].controlPoints[1],
           data[data.keys()[num]].Q[0],data[data.keys()[num]].Q[1],data[data.keys()[num]].Qres)
xlabel('x-coordinate [-]')
ylabel('y-coordinate [-]')
title('Multiple Geometries [%d bodies]' %(len(data)))
axis('equal')
grid()

