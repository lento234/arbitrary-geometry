#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Name:           main_err.py
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
windspeed = np.array([[10.], [0.]]) # x-dir and y-dir respectively
n = [10,25,50,100,200,500,1000,2000]


#==============================================================================
# Solving for various panels
#==============================================================================
err_sor = []
err_vort = []
for n_panels in n:
        
    
    #==============================================================================
    # Initialization of bodies
    #==============================================================================
               
    tower    = body(name        = 'tower',
                    shape       = ('cylinder', n_panels),
                    chord       = 1,
                    local_pitch = 0.,
                    pivot_point = [0.5, 0.])
                               
    cylinder_vort = multiBody(dict(body         = tower,
                                   location     = [0., 0.],
                                   global_pitch = 0.))
    
    cylinder_sor =  multiBody(dict(body         = tower,
                                   location     = [0., 0.],
                                   global_pitch = 0.))
                                                                          
    #==============================================================================
    # Defining Vortex and scan points
    #==============================================================================
   
    # Plotting: Mesh grid
    x,y = meshgrid(linspace(-1,1,50),linspace(-1,1,50))
    mesh = array([concatenate(x),concatenate(y)])
     
    #==============================================================================
    # Solve the panel problem      
    #==============================================================================
    
    # Source Panel
    cylinder_sor.sourcePanel_solve(freestream=windspeed)
    
    # Vortex sheet    
    cylinder_vort.vortexPanel_solve(freestream=windspeed)

    #==============================================================================
    # Calculating the Pressure coefficent
    #==============================================================================
    
    # Source Panel
    Vtot_sor = cylinder_sor.Vinduced + windspeed
    Qt_sor = Vtot_sor[0]*cylinder_sor.tangent[0] + Vtot_sor[1]*cylinder_sor.tangent[1]
    Cp_sor = 1 - (Qt_sor/sum(windspeed))**2
    
    # Vortex sheet    
    Cp_vort = 1 - ((sum(windspeed)*cylinder_vort.tangent[0] + cylinder_vort.vortex_gamma/2)/sum(windspeed))**2
    
    # Analytical Solution
    theta = rad2deg(np.linspace(np.pi,-np.pi,(n_panels)))
    Cp_analytic = 1 - 4*sin(deg2rad(theta))**2
    
    #==============================================================================
    # Plotting the pressure coefficent of various methods      
    #==============================================================================
    
    if n_panels == 100:
        figure()
        axis([-100, 100, 1.5,-3.5])
        plot(theta,Cp_analytic, label = 'Analytical')
        plot(theta,Cp_vort,'g.-', label = 'Vortex sheet')
        plot(theta,Cp_sor,'r.-', label = 'Source Panel')
        xlabel('$\theta$')
        title('Comparison: Pressure coefficent')
        ylabel('Pressure coefficient $C_p$')
        grid('on')
        legend()

    #==============================================================================
    # Calculating the error
    #==============================================================================
    err_sor.append(sum(absolute(Cp_sor-Cp_analytic))/sum(abs(Cp_analytic)))
    err_vort.append(sum(absolute(Cp_vort-Cp_analytic))/sum(abs(Cp_analytic)))


#==============================================================================
# Plotting the convergence  of 2 methods
#==============================================================================

figure()
semilogx(n,err_sor,'b.-',label='Source panel')
semilogx(n,err_vort,'g.-', label='Vortex sheet')
grid('on')
ylabel('error')
xlabel('Number of panels')
title('Convergence')
legend()
