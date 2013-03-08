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
#import inducedVelocity # to calculate the induced velocity 


#==============================================================================
# Calculate the induced velocity due to vortex strength panel
#==============================================================================
def inducedVelocity(gamma, x0, y0, x1, y1, x2, y2):
    '''
    Calculate the induced velocity due to a source term
    '''
    
    ##########!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # Transforming from global coordinates to panel coordinates
    '''
    #Figure 11.17   
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)  

    cosAlpha = (x2 - x1)/r
    sinAlpha = (y2 - y1)/r

    #x0mx1 = x0 - x1 # else, calculated twice
    #y0my1 = y0 - y1 # else, calculated twice
    
        
    # In panel coordinates (from global coordinates)    
    x  = cosAlpha*(x0-x1)  + sinAlpha*(y0-y1)
    y  = -sinAlpha*(x0-x1) + cosAlpha*(y0-y1)

    #x2mx1 = x2 - x1
    #y2my1 = y2 - y1
    
    #x1 =     
    
    #x2 = cosAlpha*(x2-x1)+ sinAlpha*(y2-y1)
    #y2 = -sinAlpha*(x2-x1) + cosAlpha*(y2-y1)
     
    # Calculating the induction at panel coordinates
   
    #Equation 11.44 and 11.45
    up = (gamma/(2*np.pi))*(np.arctan2((y-y2), (x-x2)) - np.arctan2((y-y1), (x-x1)))
    wp = -(gamma/(4*np.pi))*np.log(((x-x1)**2 + (y-y1)**2)/((x-x2)**2 + (y-y2)**2))
        
    # Transforming the velocity from panel coordinates to global coordinates
      
    #Equation 10.7
    u = cosAlpha*up - sinAlpha*wp
    w = sinAlpha*up + cosAlpha*wp
    
    return u, w
    '''
    pass

#==============================================================================
# Function to solve for the influence matrix
#==============================================================================
def solve(collocationPoint, panelStart, panelEnd, tangent):
    '''
    solve the panel method problem
    '''
    
    # Reshape data from multibody
    collocationPoint_x, collocationPoint_y,\
    panelStart_x, panelStart_y,\
    panelEnd_x, panelEnd_y,\
    tangent_x, tangent_y                  = reshapeData.toSolve(collocationPoint,
                                                                panelStart,
                                                                panelEnd,
                                                                tangent)
    
    # Creating a unit sigma source matrix for calculation
    gamma = np.ones(np.shape(collocationPoint_x)) # unit source
    
    # Calculate the induced velocity due to unit source
    u,w = inducedVelocity(gamma, collocationPoint_x, collocationPoint_y,
                          panelStart_x, panelStart_y, 
                          panelEnd_x, panelEnd_y)
    
    # Calculating the influence matrix
    A = u*tangent_x + w*tangent_y
           
    return A

#==============================================================================
# Function to evaluate the induction due to source
#==============================================================================
def evaluate(Gamma, collocationPoint, panelStart, panelEnd):
    '''
    Evaluate the source terms
    '''
    
    Sigma,\
    collocationPoint_x, collocationPoint_y,\
    panelStart_x, panelStart_y,\
    panelEnd_x, panelEnd_y              = reshapeData.toEvaluate(Gamma,
                                                                 collocationPoint,
                                                                 panelStart,
                                                                 panelEnd)
 
    # Calculating the induced velocity due to the given Source
    u,w = inducedVelocity.calc(Gamma, collocationPoint_x, collocationPoint_y,
                               panelStart_x, panelStart_y,
                               panelEnd_x, panelEnd_y)
                               
    V_vort = np.array([np.sum(u, axis=1),
                       np.sum(w, axis=1)])
                       
    return V_vort      
                                                                 
#==============================================================================
# Calculate the RHS of the problem
#==============================================================================
def RightHandSide(collocationPoint, tangent, vortex=None, freestream=[[0.],[0.]]):
    '''
    Solve the right hand side of the panel method problem
    '''
    
    if vortex is None:
        # If no point vortex needs to be calculated
        RHS = -(freestream[0]*tangent[0] + freestream[1]*tangent[1])
    else:
        # calculate the induced velocity due to vortex
        V_vort = vortex.inducedVelocity(collocationPoint)['collocationPoint']
            
        # Right-hand-side of the equation. 
        #       Adding the freestream and vortex induction before taking the norm
        RHS = -((V_vort[0] + freestream[0])*tangent[0] + (V_vort[1] + freestream[1])*tangent[1])
        
    return RHS
    
    
    
