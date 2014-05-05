#include <iostream>
#include <fstream>
#include <cmath> 
#include <assert.h> 
#include <string>
#include <sstream>
#include "FluidElement.h"
#include "Grid.h"




using namespace std;



int main()
{




 double pi = atan(1)*4;
 double length, width, psiTbc, psiBbc,psiLbc, psiRbc;
 
//INPUTS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 length = 1;
 width = 1;
 //specify BC
 psiTbc = psiBbc = psiLbc = psiRbc = 0;
 double omega = 0.9; // SOR parameter

for (int tolInc = 2; tolInc <= 4; tolInc++)
{
  //specify iteration tolerance
  double tol = pow(10,-3*tolInc);


  //grid refinement study
  int nx, ny;
  for (int n = 1; n <= 10; n++)
  {
   //generate a new grid
   nx = ny = 10*n;
   Grid grid(nx,ny,length,width);

   double L = grid.get("length");
   double W = grid.get("width");	
   double xkl, ykl, psiManu, rhsPsi, psiInitial;
   // set IC, psiManu and rhsPsi(psiManu)
   for (int l = 0 ; l <= ny - 1 ; l++)
   { 
    for (int k = 0 ; k <= nx - 1 ; k++)	
    {
     xkl     = grid.get("x",k,l);
     ykl     = grid.get("y",k,l);
     psiManu = sin(xkl*pi/L)*sin(ykl*pi/W); // this is a fake solution to psi
     grid.set("psiManu",k,l,psiManu);
	
     rhsPsi = -pi*pi*sin(xkl*pi/L)*sin(ykl*pi/W)*(1/L/L + 1/W/W);
     grid.set("rhsPsi",k,l,rhsPsi);

     psiInitial = cos(xkl*pi/L)*cos(ykl*pi/W);
     grid.set("psi",k,l,psiInitial);
    }

   }
   // set BC
   grid.set_psi_bc(psiLbc, psiRbc, psiTbc, psiBbc);
   // solve poisson eq
   grid.solve_poisson(tol,omega);




   // generate filename
   std::ostringstream os;
   os << "GridRefinementStudy_L:" << L << "_W:" << W << "_Tol:" << tol << ".txt";
   std::string filename = os.str();



   // Grid variables to save.
   double dx = grid.get("dx");
   double errorPsiMax = grid.get("errorPsiMax");
   double errorPsiAvg = grid.get("errorPsiAvg");

   //save file
   ofstream myfile;
   myfile.open (filename.c_str(),ios::app);
   //add header to file initially
   if (n == 1)
   {
    myfile << "nx " << "dx " << "dx*dx " <<  "errorMax " << "errorAbsAvg " << endl;
   }
    myfile << nx << " " << dx << " " <<  dx*dx << " "	 << errorPsiMax << " " << errorPsiAvg << endl;
   myfile.close();
//grid.save_grid_params("s.txt"	);

  }
 }

return 0;
}


