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
# Calculate the induced velocity due to source strength panel
#==============================================================================
def inducedVelocity(sigma, x, y, x1, y1, x2, y2):
    '''
    Calculate the induced velocity due to a source term
    '''
    
    # Transforming from global coordinates to panel coordinates
    
    #Figure 11.17   
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)  

    cosAlpha = (x2 - x1)/r
    sinAlpha = (y2 - y1)/r

    xmx1 = x - x1 # else, calculated twice
    ymy1 = y - y1 # else, calculated twice
    
    # In panel coordinates (from global coordinates)    
    xp  = cosAlpha*(xmx1)  + sinAlpha*(ymy1)
    yp  = -sinAlpha*(xmx1) + cosAlpha*(ymy1)    
    
    x2p = cosAlpha*(x2 - x1)+ sinAlpha*(y2 - y1)
     
    # Calculating the induction at panel coordinates
   
    # Defining preliminary variables    
    r1 = xp**2 + yp**2
    r2 = (xp - x2p)**2 + yp**2
        
    theta1 = np.arctan2(yp, xp) # theta1
    theta2 = np.arctan2(yp, (xp - x2p))
    
    #Equation 11.21 and 11.22
    up = (sigma/(4*np.pi))*np.log(r1/r2)
    wp = (sigma/(2*np.pi))*(theta2 - theta1)
        
    # Transforming the velocity from panel coordinates to global coordinates
      
    #Equation 10.7
    u = cosAlpha*up - sinAlpha*wp
    w = sinAlpha*up + cosAlpha*wp
    
    return u, w


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
    u,w = inducedVelocity(sigma, collocationPoint_x, collocationPoint_y,
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
    u,w = inducedVelocity(Sigma, collocationPoint_x, collocationPoint_y,
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
    
    
    
