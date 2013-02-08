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



# Importing Geometries
naca = file('Geometries/NACA0012.txt', 'r')
x,y  = [], []

for line in naca:
    a, b = line.strip().split('	')
    x.append(float(a))
    y.append(float(b))
naca.close()

# Dictionary of Geometries
# a, 0.5 above y=0 and b, 0.5 below y=0
Geometries = {'a': array([x,                    # x-coordinates
                          add(y[::-1],0.5)]),   # y-coordinates
              'b': array([x,
                          add(y[::-1],-0.5)])}

# Generating plotting Mesh
meshx, meshy = meshgrid(linspace(-1,2.5,50),linspace(-1,1,100))



# Calculating the data (from multipleGeometry module)
data = multipleGeometry.inducedVelocities(Geometries, meshx, meshy, Uinf=1.)


# Plotting Data
figure(1)
for name in data.geometries:
    plot(data.points[name][0],data.points[name][1],'k')
    #quiver(data.controlPoints[name][0],data.controlPoints[name][1],data.Qx[name],data.Qy[name],data.Qres[name])
#quiver(meshx,meshy,data.Qx,data.Qy)
contourf(meshx, meshy, data.Qres)
colorbar()
axis('scaled')
xlabel('x-coordinate [-]')
ylabel('y-coordinate [-]')
title('Multiple Geometries [%d bodies]' %(data.length))

raw_input("Press Enter to continue...")