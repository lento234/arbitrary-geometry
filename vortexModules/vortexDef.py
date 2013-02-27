# -*- coding: utf-8 -*-
"""
Case:           Arbitrary Geometry
Methodology:    Solving any geometry using Panel Method
Description:    Contains the multipleGeometry.py modules 
@author:        Lento Manickathan - 1544101
"""

import numpy as np
from inducedVelocity import vortexInduction

class vortex():
    def __init__(self, x_coordinates, y_coordinates, vortex_strength):
        '''
        Initialization
        '''
        
        # Storing the data
        try:
            self.data = np.asarray([np.concatenate(x_coordinates),
                                    np.concatenate(y_coordinates),
                                    np.concatenate(vortex_strength)])
        except:
            self.data = np.asarray([x_coordinates,
                                    y_coordinates,
                                    vortex_strength])
        
        # Creating pointers
        self.coordinates_pointer = self.data[0:2,:]
        self.strength_pointer    = self.data[2,:]
                
    def inducedVelocity(self, collocationPoint=None, also_onSelf=False):
        '''
        '''
        
        # Concatenating data: If meshgrid data
        #try:
        #    M = np.shape(collocationPoint)[1]
        #    collocationPoint = np.asarray([np.concatenate(collocationPoint[0]),
        #                            np.concatenate(collocationPoint[1])])
        #                            
        #except:
        #    pass
        
        # Calculate the induced velocity       
        V_ind = vortexInduction(self.data, collocationPoint=collocationPoint, also_onSelf=also_onSelf)
        
        #try:
        #    U = np.reshape(V_ind[0],(M,-1))
        #    W = np.reshape(V_ind[1],(M,-1))
        #    
        #    return U,W
        # except:
        return V_ind