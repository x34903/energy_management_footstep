Steven Crews
PhD Candidate, Cargnegie Mellon University
screws@andrew.cmu.edu
24 Feb 2020

*****************************************************************

This code accompanies the submission of the paper

"Energy Management Through Footstep Selection For Bipedal Robots"

to the 

2020 International Conference on Intelligent Robots and Systems (IROS)

and

Robot and Automation Society Letters (RA-L)

*****************************************************************

This code reproduces the plots used in Figures 4 "Curve of Capture" and "Figure 6 "Curve of Equal Energy." 
The code also creates the functions for behind equations 14 and 16. These functions are saved in the folder "symbolic_functions"

Em [J] pre-collision energy (energy-minus)
dh [m] change in footstep height 

Equation 14: capture point
x_cp = regression_create_x_cp( Em , dh ) ;

x_cp [m] at energy level Em and height dh, this gives the intersection with the curve of capture


Equation 16: equal energy point
x = regression_create_x_cp( Em , dh ) ;

x [m] at energy level Em and height dh, this gives the intersectio with the curve of equal energy

