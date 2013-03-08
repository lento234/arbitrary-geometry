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
n_panels  = 5

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
                chord       = 1,
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
cylinder_testCase = multiBody(dict(body         = tower,
                                   location     = [0., 0.],
                                   global_pitch = 0.))
                                                                      
#==============================================================================
# Defining Vortex and scan points
#==============================================================================
# Vortex Field
#vort         = vortex([-1.],[0.],[1.])

# Plotting: Mesh grid
x,y = meshgrid(linspace(-1,1,50),linspace(-1,1,50))
mesh = array([concatenate(x),concatenate(y)])
 
 
#==============================================================================
# Calculating the induced Velocity
#==============================================================================

# Solve the panelMethod
cylinder_testCase.vortexPanel_solve(evaluationPoints  = mesh,
                                    vortexPoints      = None,
                                    freestream        = windspeed)
                                

# Total Velocity: sum of source induction and windspeed
V_tot = cylinder_testCase.Vinduced + windspeed

'''
# Solve the panelMethod
windturbine.vortexPanel_solve(evaluationPoints        = mesh,
                                    vortexPoints      = None,
                                    freestream        = windspeed)
                  
# Total Velocity: sum of source induction and windspeed
V_tot = windturbine.Vinduced + windspeed
'''


Vres = (V_tot[0]**2 + V_tot[1]**2)**0.5
Vres[Vres>2.8] = nan

#==============================================================================
# TEST PLOTTING
#==============================================================================

# Plotting 
figure()
cylinder_testCase.plot('geometry')
#contourf(x,y,reshape(Vres,shape(x)))
quiver(mesh[0],mesh[1], V_tot[0], V_tot[1],Vres)
axis([-1,1,-1,1])
axis('scaled')
colorbar()

phi = linspace(-180+(180/n_panels),180+(180/n_panels),n_panels)
vs  = 2*np.sin(np.deg2rad(phi))*windspeed[0]

figure()
plot(0.5*cos(np.deg2rad(phi)),vs)
plot(cylinder_testCase.collocationPoint[0],cylinder_testCase.vortex_gamma,'.')
#plot()
#axis('equal')
axis([-1,1,0,2])
grid('on')

#==============================================================================
# Error
#==============================================================================
err = average((abs(vs-cylinder_testCase.vortex_gamma)*abs(cylinder_testCase.geometry[0,1:]-cylinder_testCase.geometry[0,:-1]))/(sum(abs(cylinder_testCase.geometry[0,1:]-cylinder_testCase.geometry[0,:-1]))))
print err

figure()
plot(array([4,10,20,30,40,50,60,70,80,90,100,200,500,1000]),
     array([0.46635927869,0.2737,0.1486,0.1010,0.0765,0.0615,0.05148,0.04425,0.038759,0.034526,0.03111,0.015633557165,0.00627131268198,0.0031]),'.-')
grid('on')