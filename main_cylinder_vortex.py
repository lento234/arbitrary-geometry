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

# Standard python scientific module [TEMPORARY]
from pylab import *
from numpy import *
ion() # Interactive on

#==============================================================================
# Initialization of the constants
#==============================================================================

# Control Parameters
windspeed = np.array([[1.], [0.]]) # x-dir and y-dir respectively
n_panels = 20
         
#==============================================================================
# Initialization of bodies
#==============================================================================
           
tower    = body(name        = 'tower',
                shape       = ('cylinder', n_panels),
                chord       = 1,
                local_pitch = 0.,
                pivot_point = [0.5, 0.])
                                 
cylinder = multiBody(dict(body         = tower,
                          location     = [0., 0.],
                          global_pitch = 0.))
                                            
#==============================================================================
# Defining Vortex and scan points
#==============================================================================

x,y = meshgrid(linspace(-1,1,100),linspace(-1,1,100))
mesh = array([concatenate(x),concatenate(y)])
 
#==============================================================================
# Solving the vortex panel problem
#==============================================================================

cylinder.vortexPanel_solve(evaluationPoints = mesh,
                           vortexPoints     = None,
                           freestream       = windspeed)


#==============================================================================
# Calculating the induced Velocity 
#==============================================================================
cylinder_Vtot = cylinder.Vinduced + windspeed
cylinder_Vres = sqrt(cylinder_Vtot[0]**2 + cylinder_Vtot[1]**2) 


#==============================================================================
# PLOTTING
#==============================================================================
figure()
cylinder.plot('geometry')
contourf(x,y,reshape(cylinder_Vres,shape(x)))
#quiver(mesh[0],mesh[1],cylinder_Vtot[0],cylinder_Vtot[1],cylinder_Vres)
colorbar()
grid('on')
axis('scaled')