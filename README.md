Description
Object oriented c++ program to solve 2D incompressible flow problems using stream-vorticity method. This program will be validated by simulating lid driven cavity flow, and comparing the results with literature. 

INSTRUCTIONS:
if you dont intend to install Visualization Toolkit (VTK 6.1.0), go straight to the "without_VTK" directory and see the readme.txt
	To access old data, see current folder.  This code is from 15/05/2014

if you want to run with VTK, you need to install VTK 6.1.0 and CMAKE *HIGHLY RECOMMENDED* 
Download http://www.vtk.org/VTK/resources/software.html 
Download http://www.cmake.org/
once installed, enter the "Current" directory  (Not this one, the one labeled "Current")
execute "./cmake"
edit the Main.cpp in the sources directory as required and make that output folder exists  (VTK will not make a new directory ect)
../make
./Main






Program Structure
The program is implemented in C++, using object orientated programming.  The goal was to make the structure modular so that in future, a variety of problems could be solved on many different style grids, in higher dimensions.

The spatial domain is discretized and stored as a grid object.  The grid class holds a fluid element object for each discrete point in the spatial domain, as well as any parameters specific to that particular grid.  The fluid element holds all variables of interest at that particular point in space, as well as its position in space.  All the variables stored at the fluid element object level are non-dimensionalized by the reference values stored in the grid object.  Since all the variables in the governing equations are stored in the fluid element, the equations being solved are non-dimensional.  

The data is output using the Visualization Toolkit from Kitware.  All the scalar variables, and the velocity fields are written to a VTK file at specified dump period.  The data can be directly opened in ParaView and visualized directly.  Additionally, a text file is written out at the end of each simulation and contains a text file with the u-velocity component along the vertical centerline, and the v-velocity component along the horizontal centerline.  The program has no runtime arguments.  All the parameters are specified in the Main.cpp file. For instruction on compiling and executing the program, refer to the readme.txt in the main project folder.
