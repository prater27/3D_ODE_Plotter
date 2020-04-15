# 3D_ODE_Plotter

3D_ODE_Plotter is tool developed in C++ using VTK. It allows creating animations of solutions for (up to) 3D Ordinary Differential Equations (also transitory ODEs), so dx/dt=f(x(t), t), by defining the vector field f(x(t), t) and giving the solution as a file text with 8 fields:
t x y z - - - -   The last four fields are dummy fields (I reused for reading the file a previous code I had where values for a quaternion giving a rotation was considered -see ). See out.txt as an example. This file contains a solution for a Lorentz ODE.

The tool can also be used for plotting non-stationary vector fields.

To compile it, you need to have installed VTK libraries and then run: cmake .   and make. To run the program then: ./ODEs_visualizator solutionFile

There is two branches. One draws the closest vectors to the current position/time solution present in the predefined grid using a different mapper, so it can be highlighted while the rest of the vector field remains hidden/semi-transparent. In the other branch, the actual tangent vector to the current position/time solution is drawn.

Some gifs of things that can be achieve with this tool are:

