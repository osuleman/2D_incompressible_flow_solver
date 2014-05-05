#include <iostream>
#include <fstream>
#include <cmath> 
#include <assert.h> 
#include <string>
#include <sstream>
#include "FluidElement.h"
#include "Grid.h"
#include "to_str.h"

using namespace std;



int main()
{

 double pi = atan(1)*4;
 double length, width, psiTbc, psiBbc,psiLbc, psiRbc;
 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 length = 1;
 width = 1;
 //specify BC
 psiTbc = psiBbc = psiLbc = psiRbc = 0;

 //specify poisson soln params
 double tol = pow(10,-10);
 double omega = 0.9;

 //grid refinement study
 int nx, ny;
 for (int n = 1; n <= 30; n++)
 {
 	nx = ny = 10*n;


	 	Grid grid(nx,ny,length,width);

 	double L = grid.get("length");
 	double W = grid.get("width");	
	double xkl, ykl, psiManu, rhs, psiInitial;




// set IC, psiManu and rhs(psiManu)

for (int l = 0 ; l <= ny - 1 ; l++)
{

	for (int k = 0 ; k <= nx - 1 ; k++)
	{

	
	
	xkl     = grid.get("x",k,l);
	ykl     = grid.get("y",k,l);
	psiManu = sin(xkl*pi/L)*sin(ykl*pi/W); // this is a fake solution to psi
	grid.set("psiManu",k,l,psiManu);
		
	rhs = -pi*pi*sin(xkl*pi/L)*sin(ykl*pi/W)*(1/L/L + 1/W/W);
	grid.set("rhs",k,l,rhs);

	psiInitial = cos(xkl*pi/L)*cos(ykl*pi/W);
	grid.set("psi",k,l,psiInitial);
	}

}

	grid.set_psi_bc(psiLbc, psiRbc, psiTbc, psiBbc);





// solve poisson eq
grid.solve_poisson(tol,omega);


if (n == 1)
{
grid.create_grid_params_header("test.txt");
}
grid.save_grid_params("test.txt");

grid.save_vtk_scaler_field("test.vti");




}













return 0;
}


