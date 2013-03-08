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

# Main Calculating module
#from panelMethod.panelMethod import panelMethod # Controls the panel method problem

# Standard python scientific module [TEMPORARY]
from pylab import *
from numpy import *
ion() # Interactive on

#==============================================================================
# Initialization of the constants
#==============================================================================

# Control Parameters
windspeed = [[1.], [0.]] # x-dir and y-dir respectively
n_panels  = 100.

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
                shape       = ('cylinder', 100),
                chord       = 0.25,
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

#==============================================================================
# Defining Vortex and scan points
#==============================================================================
# Vortex Field
vort         = vortex([-1.],[0.],[1.])

# Plotting: Mesh grid
x,y = meshgrid(linspace(-3,3,50),linspace(-3,3,50))
mesh = array([concatenate(x),concatenate(y)])
 
 
#==============================================================================
# Calculating the induced Velocity
#==============================================================================

# Solve the panelMethod
windturbine.sourcePanel_solve(evaluationPoints  = 'self',
                              vortexPoints      = None,
                              freestream        = windspeed)

# Total Velocity: sum of source induction and windspeed
V_tot = windturbine.source_Vinduced + windspeed
Vres = (V_tot[0]**2 + V_tot[1]**2)**0.5

#==============================================================================
# TEST PLOTTING
#==============================================================================

# Plotting 
figure()
windturbine.plot('geometry')
axis('scaled')
#axis([-3,3,-3,3])
#contourf(x,y,reshape(Vres,shape(x)))
quiver(windturbine.collocationPoint[0],windturbine.collocationPoint[1],V_tot[0],V_tot[1],Vres)
colorbar()