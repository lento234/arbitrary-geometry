Arbitrary Geometry
============

Similar evaluation as the problem of potential flow evaluation of single cylinder where the cylinder was defined using constant strength source panels.

Arbitrary geometry package, designed in object-oriented fashion. 

Files Include:
=============

- main.py  [main]
    call module:
        contains only scripts for initializing the problem

- multipleGeometry.py [calculation module]
	Main modules:
	    - sourceTerm:
	        calculate only the source term for the input problem. Saves the data
	        in an organized class called 'data'.
	        
	    - inducedVelocities:
	         calculates source term in necessary, and calculates the induced
	        velocities at the control points.Saves the data
	        in an organized class called 'data'.
	        
	        
    Other modules:
        - dataClass [class]: 
            creates an empty class, without attributes
            
        - Geometry [class]: 
            create a class with .x (x-axis) and .y (y-axis) data attributes
            
        - organizeData:
        
        - reshapeData2Solve:
         
        - reshapeData2Calc:
        
        
        - sor2D:

        - calc_normVec:
        
        - calc_tangVec:
        
        - calc_controlPoints:
        
        - calc_panelPoints:
        
        - panel2global
                
        - global2panel:
        
	
	
- plotMod [plotting module] ??? (not written yet)

External Modules:
- numpy

Future Work:
===========

- include calculation of induction to random field of control points, if geometries and source strengths are known. [work level = 5]
- include calculation of random vortex field vortex field. [work level = 8]

