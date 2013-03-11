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

# Vortex modules
from vortexModules.vortexDef import vortex # Contains modules related to vortex

# Standard python scientific module [TEMPORARY]
from pylab import *
from numpy import *
ion() # Interactive on

#==============================================================================
# Initialization of the constants
#==============================================================================

# Control Parameters
windspeed = np.array([[10.], [0.]]) # x-dir and y-dir respectively
n_panels = 100
         
#==============================================================================
# Initialization of bodies
#==============================================================================

# Creating bodies
airfoilA = body(name        = 'airfoilA',
                shape       = ('geometries/NACA0012.txt', False), 
                chord       = 1.,
                local_pitch = 5.,
                pivot_point = [0.25, 0.])
        
airfoilB = body(name        = 'airfoilB',
                shape       = ('geometries/NACA0012.txt', False),
                chord       = 1.,
                local_pitch = 5.,
                pivot_point = [0.25, 0.])
           
tower    = body(name        = 'tower',
                shape       = ('cylinder', n_panels),
                chord       = 1,
                local_pitch = 0.,
                pivot_point = [0.5, 0.])

# Creating multi-body
windturbine = multiBody(dict(body         = tower,
                             location     = [0.,0.],
                             global_pitch = 0.),
                             
                        dict(body         = airfoilA,
                             location     = [0.,2.],
                             global_pitch = 0.),  
                        
                        dict(body         = airfoilB,
                             location     = [0.,-2.], 
                             global_pitch = 180.))                      
                                   
cylinder = multiBody(dict(body         = airfoilA,
                          location     = [0., 0.],
                          global_pitch = 0.))
                                            
#==============================================================================
# Defining Vortex and scan points
#==============================================================================
# Vortex Field
#vort         = vortex([-1.],[0.],[1.])

# Plotting: Mesh grid
x,y = meshgrid(linspace(-1,1,100),linspace(-1,1,100))
mesh = array([concatenate(x),concatenate(y)])
x2,y2 = meshgrid(linspace(-2,2,100),linspace(-3,3,100))
mesh2 = array([concatenate(x2),concatenate(y2)])
 
#==============================================================================
# Calculating the induced Velocity
#==============================================================================

windturbine.vortexPanel_solve(evaluationPoints = mesh2,
                              vortexPoints     = None,
                              freestream       = windspeed)

cylinder.vortexPanel_solve(evaluationPoints = mesh,
                           vortexPoints     = None,
                           freestream       = windspeed)


windturbine_Vtot = windturbine.Vinduced + windspeed
windturbine_Vres = sqrt(windturbine_Vtot[0]**2 + windturbine_Vtot[1]**2) 

cylinder_Vtot = cylinder.Vinduced + windspeed
cylinder_Vres = sqrt(cylinder_Vtot[0]**2 + cylinder_Vtot[1]**2) 


#==============================================================================
# PLOTTING
#==============================================================================

figure()
cylinder.plot('geometry')
contourf(x,y,reshape(cylinder_Vres,shape(x)))
quiver(mesh[0],mesh[1],cylinder_Vtot[0],cylinder_Vtot[1],cylinder_Vres)
colorbar()
grid('on')
axis('scaled')

figure()
windturbine.plot('geometry')
contourf(x2,y2,reshape(windturbine_Vres,shape(x2)))
#quiver(mesh2[0],mesh2[1],windturbine_Vtot[0],windturbine_Vtot[1],windturbine_Vres)
colorbar()
grid('on')
axis('scaled')
