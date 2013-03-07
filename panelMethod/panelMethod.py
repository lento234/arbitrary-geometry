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

#==============================================================================
# Main panel method evaluation module
#==============================================================================
def sourcePanel(collocationPoint, panelStart, panelEnd, normal, evaluationPoint, vortexParticles=None, freestream = [[0.],[0.]]):
    '''
    Solves the potential flow problem using source terms
    '''
    
    # Calculate RHS - Neumann B.C
    RHS = source2D.RightHandSide(collocationPoint, 
                                 normal, 
                                 vortexParticles, 
                                 freestream)
    
    # Calculate the influence matrix
    A   = source2D.solve(collocationPoint, 
                         panelStart,
                         panelEnd,
                         normal)
    
    # Solve the panel Method. Equation: Ax = RHS. Solve for x
    Sigma = np.linalg.solve(A,RHS)

    # To show no transpiration    
    #if evaluationPoint is 'self':
    #    evaluationPoint = collocationPoint
        
    # Calculate induced velocity on body
    V_sorc = source2D.evaluate(Sigma,
                               evaluationPoint,
                               panelStart,
                               panelEnd)
    
    return RHS, A, Sigma, V_sorc
    
    
def vortexPanel(collocationpoint, panelStart, panelEnd, tangent, evaluationPoint, vortexParticles=None, freestream= [[0.],[0.]]):
    '''
    Solves the potential flow problem using vortex panels
    '''
    
    pass
    
    
    