# MATLAB-Lid-Driven-Cavity-Flow
SImulation of translating lid flow in a square cavity using unsteady 2D Navier-Stokes equations  
  
This simulation simulates a translating lid flow in a square cavity. To do this the nondimensional incompressible unsteady 2D Navier-Stokes equations are used. The method used is the finite volume method on a staggered cartesian grid with square uniform cells with the central difference scheme applied. The advantage of using a staggered grid is an improved accuracy and stability of the algorithm
avoiding odd-even decoupling and non-physical solutions while the disadvantage is to keep track of the different grids used.  
  
Note that the simulation might run for some time (1-2min).
