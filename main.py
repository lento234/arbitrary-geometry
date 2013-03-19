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
n_panels = 100


#==============================================================================
# Initialization of bodies
#==============================================================================

# Creating bodies
airfoil = body(name         = 'airfoil',
               shape        = ('geometries/NACA0012.txt', False), 
               chord        = 1.,
               local_pitch  = 0.,
               pivot_point  = [0.25, 0.])
           
tower   = body(name         = 'tower',
               shape        = ('cylinder', n_panels),
               chord        = 1.,#sqrt(1.+1.),
               local_pitch  = 0.,
               pivot_point  = [0.5, 0.])
               
box     = body(name         = 'block',
               shape        = ('block', n_panels),
               chord        = 1.,
               local_pitch  = 0.,
               pivot_point  = [0., 0.])
               

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
                          global_pitch = -360./(2*n_panels))) 

#block   = multiBody(dict(body         = box,
#                         location     = [0., 0.],
#                         global_pitch = 0.))
                                                             
#==============================================================================
# Defining Vortex and scan points
#==============================================================================
# Vortex Field
#vort         = vortex([-1.],[0.],[1.])

# Plotting: Mesh grid
x,y = meshgrid(linspace(-1,1,200),linspace(-1,1,200))
mesh = array([concatenate(x),concatenate(y)])
 
V_external           = array([repeat(windspeed[0],shape(cylinder.collocationPoint)[1]),
                              repeat(windspeed[1],shape(cylinder.collocationPoint)[1])])

#==============================================================================
# Calculating the induced Velocity
#==============================================================================

# Source Panel                              
cylinder.sourcePanel_solve(Vinduced         = V_external,
                           evaluationPoint  = mesh)                                               

Vtot_source = cylinder.source_V + windspeed
Vres_source = (Vtot_source[0]**2 + Vtot_source[1]**2)**0.5


# Vortex panel
cylinder.vortexPanel_solve(Vinduced         = V_external,
                           evaluationPoint  = mesh)

Vtot_vortex = cylinder.vortex_V + windspeed
Vres_vortex = (Vtot_vortex[0]**2 + Vtot_vortex[1]**2)**0.5                                               
                                              

# Source-Vortex Panel         
cylinder.sourceVortexPanel_solve(Vinduced         = V_external,
                                 evaluationPoint  = mesh)  

Vtot_sourceVortex = cylinder.sourceVortex_V + windspeed
Vres_sourceVortex = (Vtot_sourceVortex[0]**2 + Vtot_sourceVortex[1]**2)**0.5

#==============================================================================
# Calculating the analytical solution
#==============================================================================

# Defining the r and theta
r = (mesh[0]**2 + mesh[1]**2)**0.5
theta = np.arctan2(mesh[1],mesh[0])

Vr      = windspeed[0]*(1 - 0.5**2/(r**2))*np.cos(theta)
Vtheta  = -windspeed[0]*(1 + 0.5**2/(r**2))*np.sin(theta) 

Vres_analytical = (Vr**2 + Vtheta**2)**0.5
Vres_analytical[r<0.5] = nan
#Vx      = np.cos(theta)*Vr - r*np.sin(theta)*Vtheta
#Vy      = np.sin(theta)*Vr + r*np.cos(theta)*Vtheta

# Calculating the error
diff_source         = abs((Vres_analytical - Vres_source))#/Vres_analytical)
diff_vortex         = abs((Vres_analytical - Vres_vortex))#/Vres_analytical)
diff_sourceVortex   = abs((Vres_analytical - Vres_sourceVortex))#/Vres_analytical)

#==============================================================================
# Plots
#==============================================================================

figure()
cylinder.plot('geometry')
contourf(x,y,reshape(Vres_analytical,shape(x)))
#quiver(mesh[0],mesh[1],Vtot[0],Vtot[1],Vres)
#plot(x,y,'k.')
axis('scaled')
grid('on')
colorbar()
title('Analytical Solution')
xlabel('x')
ylabel('y')

figure()
cylinder.plot('geometry')
contourf(x,y,reshape(Vres_source,shape(x)))
axis('scaled')
grid('on')
colorbar()
title('Source')
ylabel('y')
xlabel('x')

figure()
cylinder.plot('geometry')
contourf(x,y,reshape(Vres_vortex,shape(x)))
axis('scaled')
grid('on')
colorbar()
title('Vortex')
ylabel('y')
xlabel('x')

figure()
cylinder.plot('geometry')
contourf(x,y,reshape(Vres_sourceVortex,shape(x)))
axis('scaled')
grid('on')
colorbar()
title('Source-Vortex')
ylabel('y')
xlabel('x')

figure()
cylinder.plot('geometry')
contourf(x,y,reshape(diff_source,shape(x)))
axis('scaled')
grid('on')
colorbar()
title('Error - Source')
ylabel('y')
xlabel('x')

figure()
cylinder.plot('geometry')
contourf(x,y,reshape(diff_vortex,shape(x)))
axis('scaled')
grid('on')
colorbar()
title('Error - Vortex')
ylabel('y')
xlabel('x')

figure()
cylinder.plot('geometry')
contourf(x,y,reshape(diff_sourceVortex,shape(x)))
axis('scaled')
grid('on')
colorbar()
title('Error - Source Vortex')
ylabel('y')
xlabel('x')
#quiver(mesh[0],mesh[1],Vx,Vy)

#xlabel('$\theta$')
#title('Comparison: Pressure coefficent')
#ylabel('Pressure coefficient $C_p$')
