#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Name:           main.py
Description:    Main command module. 
Author:         Lento Manickathan - 1544101 - lento.manickathan@gmail.com
"""

#==============================================================================
# Importing Modules
#==============================================================================

# Body modules
from bodyModules.body import body # Stores all the body parameters
from bodyModules.multiBody import multiBody # Stores multiple geometries
from bodyModules.getCoordinates import getCoordinates # Extracts coordinates
from bodyModules.generate2DCylinder import generate2DCylinder # Create Circle

# Vortex modules
from vortexModules.vortexDef import vortex # Contains modules related to vortex

# Main Calculating module
from panelMethod.panelMethod import panelMethod # Controls the panel method problem

# Standard python scientific module [TEMPORARY]
from pylab import *
from numpy import *
ion() # Interactive on

#==============================================================================
# Initialization of the constants
#==============================================================================

# Control Parameters
windspeed = [1., 0.] # x-dir and y-dir respectively
n_panels  = 100.

#==============================================================================
# Initialization of bodies and vortex field
#==============================================================================

# Extracting the coordinates from files
NACA0012_coor = getCoordinates(filename  = 'Geometries/NACA0012.txt',
                               clockwise = False) # extract data from file, 
                                                  # coordinates are given anti-clockwise       

cylinder_coor = generate2DCylinder(n = n_panels) # outputs cylinder coordinates

# Creating bodies
airfoilA = body(name        = 'airfoilA',
                coordinates = NACA0012_coor, 
                chord       = 1.,
                local_pitch = 10.,
                pivot_point = [0.25, 0.])
                
airfoilB = body(name        = 'airfoilB',
                coordinates = NACA0012_coor,
                chord       = 1.,
                local_pitch = 5.,
                pivot_point = [0.25, 0.])
                
tower    = body(name        = 'tower',
                coordinates = cylinder_coor,
                chord       = 0.5,
                local_pitch = 0.,
                pivot_point = [0.5, 0.])
                               
'''
# Creating multi-body
windturbine = multiBody(dict(body         = airfoilA,
                             location     = [0.,2.],
                             global_pitch = 0.),  
                        
                        dict(body         = airfoilB,
                             location     = [0.,-2.], 
                             global_pitch = 180.),

                        dict(body         = tower,
                             location     = [0.,0.],
                             global_pitch = 0.))



# Vortex Field
vortx, vorty = np.meshgrid(np.linspace(-2,-1,10),np.linspace(-1.,1.,10))
vortGamma    = np.random.rand(np.shape(vortx)[0],np.shape(vortx)[1])
vort         = vortex(vortx,vorty,vortGamma) 

# Plotting: Mesh grid
#x,y = meshgrid(linspace(-1.5,2.,100),linspace(-1.,1.,100))
#mesh = array([concatenate(x),concatenate(y)])
 
#==============================================================================
# Solves the Multi-Body problem         
#==============================================================================

# Solve the panelMethod
V_tot = panelMethod(bodies         = windturbine, 
                    meshField_coor = vort.coordinates_pointer, 
                    vortexField    = vort, 
                    Freestream     = windspeed)
    
#Vres = (V_tot[0]**2 + V_tot[1]**2)**0.5
'''
#==============================================================================
# TEST PLOTTING
#==============================================================================

# Plotting    
#vort.coordinates_pointer = vort.coordinates_pointer + V_tot*dt
#ax1.set_xdata(vort.coordinates_pointer[0])
#ax1.set_ydata(vort.coordinates_pointer[1])
#draw()    
#axis('scaled')
#pause(1/30.)

#quiver(meshField_coor[0], meshField_coor[1],V_tot[0],V_tot[1])
#plot(windturbine.geometry[0],windturbine.geometry[1],'b.')
#plot(vortx,vorty,'go')
#axis('scaled')
#axis([-2,5,-2,2])

#figure()
#quiver(windturbine.collocationPoint[0],windturbine.collocationPoint[1],
#       V_tot[0],V_tot[1],Vres)
#axis([-0.5,1,-0.5,0.5])
#axis('scaled')


#figure()
#quiver(windturbine.collocationPoint[0],windturbine.collocationPoint[1],
#     windturbine.normal[0],windturbine.normal[1])
#axis('scaled')
#V_tot[V_tot>20.] = nan
#Vres[Vres>20.] = nan
#figure()
#plot(vort.coordinates_pointer[0],vort.coordinates_pointer[1],'o')
#quiver(vort.coordinates_pointer[0],vort.coordinates_pointer[1], V_tot[0], V_tot[1],Vres)
##contourf(x,y,Vtot)
#axis('scaled')
#colorbar()

'''
pl.figure()
pl.plot(a.geometry[0],a.geometry[1],'b.-')
pl.plot(a.collocationPoint[0],a.collocationPoint[1],'go')
pl.plot(vort.coordinates_pointer[0],vort.coordinates_pointer[1],'ro')
pl.axis('scaled')
#pl.axis([-2,2,-2,2])


V = vort.inducedVelocity(a.collocationPoint)
pl.quiver(a.collocationPoint[0],a.collocationPoint[1],V['collocationPoint'][0],V['collocationPoint'][1])

pl.show()
#x,y = meshgrid(linspace(-3,3,50),linspace(-3,3,50))

'''


#x = concatenate(x)
#y = concatenate(y)

#V_mesh = vort.inducedVelocity((x,y))
#V_mesh = reshape(average(V_mesh['collocationPoint'],axis=0),(50,50))
#x,y = reshape(x,(50,50)), reshape(y,(50,50))

#quiver(a.collocationPoint[0],a.collocationPoint[1],V['collocationPoint'][0],V['collocationPoint'][1])


#figure()
#plot(a.geometry[0],a.geometry[1],'k.-')
#plot(vort.coordinates_pointer[0],vort.coordinates_pointer[1],'go')
#quiver(vortx,vorty,V['vortex'][0],V['vortex'][1],average(V['vortex'],axis=0))
#contourf(x,y,V_mesh)
#colorbar()



# Solving the data
#figure()
#plot(airfoilA.geometry_normalized[0],airfoilA.geometry_normalized[1])
#figure()
#plot(wt.geometry[0],wt.geometry[1],'.')
#plot(wt.collocationPoint[0],wt.collocationPoint[1],'r.')
#axis('scaled')
#grid('on')

#figure()
#plot(airfoilA.geometry[0],airfoilA.geometry[1])
#plot(0.,0,'o')
#axis('scaled')
#grid('on')