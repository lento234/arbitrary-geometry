# -*- coding: utf-8 -*-
"""

"""
#==============================================================================
# Importing Modules
#==============================================================================

# Standard scientific module
import numpy as np

#==============================================================================
# Restructure the data array to solve for source term
#==============================================================================
def toSolve(collocationPoint, panelStart, panelEnd, normal):
    '''
    Reshape the data and concatenate to one matrix, tile to right format
    to solve for the source term. The end shape of every matrix will be
    (N x M), where N is number of points and M is number of panels.
    The control point attributes vary by rows, whereas the panel attributes
    vary by column, therby enabling the evaluation of the influence of
    every panel to every control points.
    
    Parameters
    ----------
    self - the output from ._concatenateData().
    
    Returns
    -------
    Reshaped variables. Shape: (N by M). N = Number of control Points,
                            M = Number of panels. The control point 
                            attributes vary by rows, whereas the panel 
                            attributes vary by column, therby enabling the 
                            evaluation of the influence of every panel to 
                            every control points. 
    
    '''
    
    # The end shape of the data.
    N = np.shape(collocationPoint)[1] # total number of control points
    M = np.shape(panelStart)[1] # total number of panel points
  
    # Control points vary row-wise
    collocationPoint_x = np.tile(np.array([collocationPoint[0]]).transpose(), [1,M])
    collocationPoint_y = np.tile(np.array([collocationPoint[1]]).transpose(), [1,M])
    
    # Panel points vary column-wise
    panelStart_x = np.tile(panelStart[0], [N,1])
    panelStart_y = np.tile(panelStart[1], [N,1])
    
    panelEnd_x = np.tile(panelEnd[0], [N,1])
    panelEnd_y = np.tile(panelEnd[1], [N,1])
    
    # UnitVectors are attributes of THE PANEL AT THE CONTROL POINT. 
    #       ergo, varies row-wise aswell, same as control points
    normal_x = np.tile(np.array([normal[0]]).transpose(),[1,M])
    normal_y = np.tile(np.array([normal[1]]).transpose(),[1,M])
   
    return collocationPoint_x, collocationPoint_y, panelStart_x, panelStart_y,\
           panelEnd_x, panelEnd_y, normal_x, normal_y#, tangent_x, tangent_y
  
#==============================================================================
# Restructure the data array to evalute the induction due to source term         
#==============================================================================
def toEvaluate(Sigma, collocationPoint, panelStart, panelEnd):
    '''
    Reshape the data for evaluate the source 
    '''
    
    # The end shape of the data.
    N = np.shape(collocationPoint)[1] # total number of control points
    M = np.shape(panelStart)[1] # total number of panel points
    
    # Control points vary row-wise
    collocationPoint_x = np.tile(np.array([collocationPoint[0]]).transpose(), [1,M])
    collocationPoint_y = np.tile(np.array([collocationPoint[1]]).transpose(), [1,M])
    
    # Panel points vary column-wise
    panelStart_x = np.tile(panelStart[0], [N,1])
    panelStart_y = np.tile(panelStart[1], [N,1])
    
    panelEnd_x = np.tile(panelEnd[0], [N,1])
    panelEnd_y = np.tile(panelEnd[1], [N,1])
  
    Sigma = np.tile(Sigma, [N,1]) # Source term vary column-wise
    
    return Sigma, collocationPoint_x, collocationPoint_y,\
           panelStart_x, panelStart_y,  panelEnd_x, panelEnd_y
    
    
           
           