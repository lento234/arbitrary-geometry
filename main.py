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
windspeed = np.array([[10.], [0.]]) # x-dir and y-dir respectively
n = range(10,500,50)

err_sor = []
err_vort = []
for n_panels in n:
        
    
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
                                       
    cylinder_vort = multiBody(dict(body         = tower,
                                   location     = [0., 0.],
                                   global_pitch = 0.))
    
    cylinder_sor =  multiBody(dict(body         = tower,
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
    
    
    cylinder_sor.sourcePanel_solve(freestream=windspeed)
    
    Vtot_sor = cylinder_sor.Vinduced + windspeed
    Qt_sor = Vtot_sor[0]*cylinder_sor.tangent[0] + Vtot_sor[1]*cylinder_sor.tangent[1]
    Cp_sor = 1 - (Qt_sor/sum(windspeed))**2
    
    
    cylinder_vort.vortexPanel_solve(freestream=windspeed)
    Cp_vort = 1 - ((sum(windspeed)*cylinder_vort.tangent[0] + cylinder_vort.vortex_gamma/2)/sum(windspeed))**2
    
    theta = rad2deg(np.linspace(np.pi,-np.pi,(n_panels)))
    Cp_analytic = 1 - 4*sin(deg2rad(theta))**2
    

    err_sor.append(sum(absolute(Cp_sor-Cp_analytic))/sum(abs(Cp_analytic)))
    err_vort.append(sum(absolute(Cp_vort-Cp_analytic))/sum(abs(Cp_analytic)))


figure()
semilogx(n,err_sor,'b.-')
semilogx(n,err_vort,'g.-')
grid('on')
