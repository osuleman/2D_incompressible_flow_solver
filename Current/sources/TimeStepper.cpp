#include <iostream>
#include <fstream>
#include <cmath> 
#include <assert.h> 
#include <string>
#include <sstream>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>

#include "FluidElement.h"
#include "Grid.h"
#include "TimeStepper.h"

using namespace std;

TimeStepper::TimeStepper(double tiSpec, double tfSpec, int nxSpec,int nySpec, double lSpec, double wSpec, double lRefSpec, double uRefSpec)
{

ti = tiSpec;
tf = tfSpec;
grid.set_grid(nxSpec,nySpec,lSpec,wSpec,lRefSpec,uRefSpec);
grid.t = ti;

}


void TimeStepper::Run(double cflSpec, int dumpPeriodSpec)
{

 if (grid.dx >= grid.dy)
 {
  dt = cflSpec*grid.dx/grid.uRef;
 }
else
 {
  dt = cflSpec*grid.dy/grid.uRef;
 }

Grid newGrid = grid;



}





