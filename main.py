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
#from vortexModules.vortexDef import vortex # Contains modules related to vortex

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
windspeed = np.array([[1.], [0.]]) # x-dir and y-dir respectively
n_panels = 18
        

#==============================================================================
# Initialization of bodies
#==============================================================================

# Creating bodies
airfoilA = body(name        = 'airfoilA',
                shape       = ('geometries/NACA0012.txt', False), 
                chord       = 1.,
                local_pitch = 0.,
                pivot_point = [0.25, 0.])
        
airfoilB = body(name        = 'airfoilB',
                shape       = ('geometries/NACA0012.txt', False),
                chord       = 1.,
                local_pitch = 0.,
                pivot_point = [0.25, 0.])
           
tower    = body(name        = 'tower',
                shape       = ('cylinder', n_panels),
                chord       = 1,#sqrt(1.+1.),
                local_pitch = 0.,
                pivot_point = [0.5, 0.])
'''
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
'''                      
                                   
cylinder = multiBody(dict(body         = tower,
                          location     = [0., 0.],
                          global_pitch = 0.))#-360./(2*n_panels))) 
                                                             
#==============================================================================
# Defining Vortex and scan points
#==============================================================================
# Vortex Field
#vort         = vortex([-1.],[0.],[1.])

# Plotting: Mesh grid
x,y = meshgrid(linspace(-1,1,200),linspace(-1,1,200))
mesh = array([concatenate(x),concatenate(y)])
 
#==============================================================================
# Calculating the induced Velocity
#==============================================================================

V_external           = array([repeat(windspeed[0],shape(cylinder.collocationPoint)[1]),
                              repeat(windspeed[1],shape(cylinder.collocationPoint)[1])])
                              
#cylinder.vortexPanel_solve(Vinduced         = V_external,
#                           evaluationPoint  = mesh)
                                               
#cylinder.sourcePanel_solve(Vinduced         = V_external,
#                           evaluationPoint  = mesh)                                               
                                               
cylinder.sourceVortexPanel_solve(Vinduced         = V_external,
                                 evaluationPoint  = mesh)
                           
Vtot = cylinder.sourceVortex_V + windspeed
Vres = (Vtot[0]**2 + Vtot[1]**2)**0.5

#==============================================================================
# Plots
#==============================================================================
figure()
cylinder.plot('geometry')
contourf(x,y,reshape(Vres,shape(x)))
#contourf(x,y,reshape(Vtot[0],shape(x)))
#quiver(mesh[0],mesh[1],Vtot[0],Vtot[1],Vres)
#plot(x,y,'k.')
axis('scaled')
grid('on')
colorbar()
#xlabel('$\theta$')
#title('Comparison: Pressure coefficent')
#ylabel('Pressure coefficient $C_p$')
