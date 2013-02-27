# -*- coding: utf-8 -*-
"""
Name:       Bodies. [Body and MultiBody modules]
Description:Contains all the data of body and its function
Author:     Lento Manickathan - 1544101
"""

#==============================================================================
# Importing modules
#==============================================================================

# Generic Scientific module
import numpy as np
# Custom - panel method modules
import reshapeData # to reshape the given datas
import inducedVelocity # to calculate the induced velocity 


#==============================================================================
# Function to solve for the influence matrix
#==============================================================================
def solve(collocationPoint, panelStart, panelEnd, normal):
    '''
    solve the panel method problem
    '''
    
    # Reshape data from multibody
    collocationPoint_x, collocationPoint_y,\
    panelStart_x, panelStart_y,\
    panelEnd_x, panelEnd_y,\
    normal_x, normal_y                  = reshapeData.toSolve(collocationPoint,
                                                              panelStart,
                                                              panelEnd,
                                                              normal)
    
    # Creating a unit sigma source matrix for calculation
    sigma = np.ones(np.shape(collocationPoint_x)) # unit source
    
    # Calculate the induced velocity due to unit source
    u,w = inducedVelocity.calc(sigma, collocationPoint_x, collocationPoint_y,
                               panelStart_x, panelStart_y, 
                               panelEnd_x, panelEnd_y)
    
    # Calculating the influence matrix
    A = u*normal_x + w*normal_y
           
    return A

#==============================================================================
# Function to evaluate the induction due to source
#==============================================================================
def evaluate(Sigma, collocationPoint, panelStart, panelEnd):
    '''
    Evaluate the source terms
    '''
    
    Sigma,\
    collocationPoint_x, collocationPoint_y,\
    panelStart_x, panelStart_y,\
    panelEnd_x, panelEnd_y              = reshapeData.toEvaluate(Sigma,
                                                                 collocationPoint,
                                                                 panelStart,
                                                                 panelEnd)
 
    # Calculating the induced velocity due to the given Source
    u,w = inducedVelocity.calc(Sigma, collocationPoint_x, collocationPoint_y,
                               panelStart_x, panelStart_y,
                               panelEnd_x, panelEnd_y)
                               
    V_sorc = np.array([np.sum(u, axis=1),
                       np.sum(w, axis=1)])
                       
    return V_sorc      
                                                                 
#==============================================================================
# Calculate the RHS of the problem
#==============================================================================
def RightHandSide(collocationPoint, normal, vortex = None, Freestream = [0., 0.]):
    '''
    Solve the right hand side of the panel method problem
    '''
    
    if vortex is None:
        # If no point vortex needs to be calculated
        RHS = -(Freestream[0]*normal[0] + Freestream[1]*normal[1])
    else:
        # calculate the induced velocity due to vortex
        V_vort = vortex.inducedVelocity(collocationPoint)['collocationPoint']
            
        # Right-hand-side of the equation. 
        #       Adding the freestream and vortex induction before taking the norm
        RHS = -((V_vort[0] + Freestream[0])*normal[0] + (V_vort[1] + Freestream[1])*normal[1])
        
    return RHS
    
    
    
