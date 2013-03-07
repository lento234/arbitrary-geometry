# -*- coding: utf-8 -*-
"""
Name:       Bodies. [Body and MultiBody modules]
Description:Contains all the data of body and its function
Author:     Lento Manickathan - 1544101
"""

#==============================================================================
# Importing necessary modules 
#==============================================================================

# Constant-strength source potential flow module
import source2D
# Standard python-scientific module
import numpy as np
#import pylab as py


#==============================================================================
# Main panel method evaluation module
#==============================================================================
def panelMethod(bodies, meshField_coor='self', vortexField=None, freestream = [0.,0.]):
    '''
    Solves the potential flow problem using source terms
    '''
    
    # Calculate RHS - Neumann B.C
    RHS = source2D.RightHandSide(bodies.collocationPoint, 
                                 bodies.normal, 
                                 vortexField, 
                                 Freestream)
    
    # Calculate the influence matrix
    A   = source2D.solve(bodies.collocationPoint, 
                         bodies.panelStart,
                         bodies.panelEnd,
                         bodies.normal)
    
    # Solve the panel Method. Equation: Ax = RHS. Solve for x
    Sigma = np.linalg.solve(A,RHS)

    # To show no transpiration    
    if meshField_coor is 'self':
        meshField_coor = bodies.collocationPoint
        
    # Calculate induced velocity on body
    V_sorc = source2D.evaluate(Sigma,
                               meshField_coor,
                               bodies.panelStart,
                               bodies.panelEnd)
    
    if vortexField is None:
        V_tot = V_sorc + np.array([Freestream]).transpose()
    else:
        # Calculate the induced velocity due to vortex field
        V_vort = vortexField.inducedVelocity(meshField_coor)
        # Total velocity
        V_tot = V_sorc + V_vort['collocationPoint'] + np.array([Freestream]).transpose()
    
    return V_tot
    
    
    
    
    