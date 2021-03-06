#include <iostream>
#include <fstream>
#include <cmath> 
#include <assert.h> 
#include <string>
#include <sstream>
#include "FluidElement.h"
#include "Grid.h"
#include "GridRefinementStudy.h"


using namespace std;



int main()
{

 double pi = atan(1)*4;

 // solving on 3 grids
 int    gridS[3]  = {33,65,129};
 int    tolExp[3] = {-7,-8,-9};
 double omegaS[3]  = {1.82,1.91,1.95}; 


 //solving for 2 different reynolds numbers
 double reynoldS[2] = {100,1000};

 //lid velocity
 double vLbc = 0;
 double vRbc = 0;
 double uTbc = 1;
 double uBbc = 0;
 double uRef = 1;

 //dims
 double l    = 1;
 double w    = 1;
 double lRef = l;


 double ti = 0;
 double tf = 0.05;
 double cfl = 0.2;
 double dumpPeriod = 500;
 std::string filenameSpec = "Results/Sim3/"; // main data storage directory
 std::string solver = "sor"; 
 std::string filename;
 std::ostringstream os;

 double L, W, xkl, ykl, psiManu, rhs, psiInitial,zetaInitial,omega,tol,reynolds;
 int N, iterations;
 Grid grid;

 for (int r = 1; r <=1; r++)
 {
  reynolds = reynoldS[r];

  for (int i = 2; i <= 2; i++)
  {

   N     = gridS[i];
   tol   = pow(10,tolExp[i]);
   omega = omegaS[i];

  // place in correct folder
  os.clear();
  os.str("");
  os << filenameSpec << "Re" << reynolds << "/N" << N << "/" ;
  filename = os.str();

   grid.set_grid(N,N,l,w,lRef,uRef);
   grid.set_boundary_velocities(vLbc,vRbc,uTbc,uBbc);

   L  = grid.get("length");
   W  = grid.get("width");  

   // set IC, psiManu and rhs(psiManu)
   for (int l = 0 ; l <= N - 1 ; l++)
   { 
    for (int k = 0 ; k <= N - 1 ; k++)	
    {
     xkl     = grid.get("x",k,l);
     ykl     = grid.get("y",k,l);

     zetaInitial = 0;
     grid.set("zeta",k,l,zetaInitial);

     psiInitial = 0;
     grid.set("psi",k,l,psiInitial);
    } 
   }// end of iteration through nodes

   // set BC
   grid.set_psi_bc("diri",0.0,0.0,0.0,0.0);
   grid.set_zeta_bc();
   grid.solve(ti,tf,reynolds,cfl,tol,omega,dumpPeriod,solver,filename);
 


 
  }// end of grid refinement loop
 }// end of Re loop


}
/*
//compare runtime
 double pi = atan(1)*4;

 // solving on 3 grids
 int    gridS[3]  = {33,65,129};
 int    tolExp[3] = {-7,-8,-9};
 double omegaS[3]  = {1.82,1.91,1.95}; 


 //solving for 2 different reynolds numbers
 double reynoldS[2] = {100,1000};

 //lid velocity
 double vLbc = 0;
 double vRbc = 0;
 double uTbc = 1;
 double uBbc = 0;
 double uRef = 1;

 //dims
 double l    = 1;
 double w    = 1;
 double lRef = l;


 double ti = 0;
 double tf = 0.05;
 double cfl = 1;
 double dumpPeriod = 100;
 std::string filenameSpec = "Results/Sim3/"; // main data storage directory
 std::string solver = "adi"; 
 std::string filename;
 std::ostringstream os;

 double L, W, xkl, ykl, psiManu, rhs, psiInitial,zetaInitial,omega,tol,reynolds;
 int N, iterations;
 Grid grid;

 for (int r = 0; r <=0; r++)
 {
  reynolds = reynoldS[r];

  for (int i = 2; i <= 2; i++)
  {

   N     = gridS[i];
   tol   = pow(10,tolExp[i]);
   omega = omegaS[i];

  // place in correct folder
  os.clear();
  os.str("");
  os << filenameSpec << "Re" << reynolds << "/N" << N << "/" ;
  filename = os.str();

   grid.set_grid(N,N,l,w,lRef,uRef);
   grid.set_boundary_velocities(vLbc,vRbc,uTbc,uBbc);

   L  = grid.get("length");
   W  = grid.get("width");  

   // set IC, psiManu and rhs(psiManu)
   for (int l = 0 ; l <= N - 1 ; l++)
   { 
    for (int k = 0 ; k <= N - 1 ; k++)	
    {
     xkl     = grid.get("x",k,l);
     ykl     = grid.get("y",k,l);

     zetaInitial = 0;
     grid.set("zeta",k,l,zetaInitial);

     psiInitial = 0;
     grid.set("psi",k,l,psiInitial);
    } 
   }// end of iteration through nodes

   // set BC
   grid.set_psi_bc("diri",0.0,0.0,0.0,0.0);
   grid.set_zeta_bc();
   grid.solve(ti,tf,reynolds,cfl,tol,omega,dumpPeriod,solver,filename);
 


 
  }// end of grid refinement loop
 }// end of Re loop


}

/*
// Optimal SOR parameter
 
  ofstream myfile;
  std::string filenameSpec = "Results/GRS/SOR_OPT";
  myfile.open (filenameSpec.c_str(),ios::app);
  myfile << "omega " << "iterations" << endl;
  myfile.close();

 int gridSize[3] = {33,65,129};
 double pi = atan(1)*4;
 double L,W,xkl, ykl, psiManu, rhs, psiInitial,omega;
 int N, iterations;
 double tol = pow(10,-9);
 double omegaI   = 0.5;
 double dOmega   = 0.01;
 int    nOmega   = 150;  
 Grid   grid;
 for (int n = 0; n <= 2; n++)
  {
   N = gridSize[n];
   for (int i = 0; i <= nOmega-1; i ++)
   {
    omega = omegaI + i*dOmega;
    grid.set_grid(N,N,1.0,1.0,1.0,1.0);
  
    L  = grid.get("length");
    W  = grid.get("width");

    // set IC, psiManu and rhs(psiManu)
    for (int l = 0 ; l <= N - 1 ; l++)
    { 
     for (int k = 0 ; k <= N - 1 ; k++)	
     {
      xkl     = grid.get("x",k,l);
      ykl     = grid.get("y",k,l);
     
      psiManu = sin(xkl*pi/L)*sin(ykl*pi/W); // this is a fake solution to psi
      grid.set("psiManu",k,l,psiManu);

      rhs = -pi*pi*sin(xkl*pi/L)*sin(ykl*pi/W)*(1/L/L + 1/W/W);
      grid.set("zeta",k,l,-rhs);

      psiInitial = cos(xkl*pi/L)*cos(ykl*pi/W);
      grid.set("psi",k,l,psiInitial);
     }
   }// end of iteration through nodes
   // set BC
   grid.set_psi_bc("diri",0, 0, 0, 0);
   // solve poisson eq
   grid.solve_poisson("sor",tol,omega);
   grid.calc_psi_error();
 
   iterations = grid.get("iterations");
   myfile.open (filenameSpec.c_str(),ios::app);
   myfile << omega << " " << iterations << endl;
   myfile.close();

  }

  myfile.open (filenameSpec.c_str(),ios::app);
  myfile << N << endl;
  myfile.close();




 }

}
*/

/*
// Grid Refinement Study

int nMin = 10;
int nGrids = 15;
int dn  = 10;
GridRefinementStudy grs(nMin,nGrids,dn);
grs.test_derivative1();
grs.test_derivative2();



// testing different SOR Parameters, omega

int    nOmega = 1;
double omegaI = 1;
double dOmega = 0.3;
double omega  = 0;

int    nTol = 7;
double expI = 6;
double dExp = 1;
double exp  = 0;
double tol  = 0;





for (int i = 0; i <= nOmega-1; i ++)
{
	for (int n = 0; n <= nTol-1; n++)
	{ 
	 exp = expI + n*dExp;
	 tol = pow(10,-exp);
         omega = omegaI + i*dOmega;
	 cout << "tol " << tol << " omega " << omega << endl;
	 grs.test_poisson_sor(tol,omega);
	 grs.test_poisson_adi(tol);
	}


}
*/
/*
//grid params
int nx        = 10;
int ny        = 10;
double length = 1;
double width  = 1;
double lRef   = length;

double ti = 0;
double tf = 3;
double reynolds = 1000;
double cfl = 0.1;
double tol = pow(10,-6);
double omega = 1;
double dumpPeriod = 1;
std::string filename = "Results/Trends/adi";

double vLbc, vRbc, uTbc, uBbc;

//lid velocity
vLbc = 0;
vRbc = 0;
uTbc = 1;
uBbc = 0;
double uRef   = 1;


Grid grid(nx,ny,length,width,lRef,uRef);
grid.set_boundary_velocities(vLbc,vRbc,uTbc,uBbc);

   double L = grid.get("length");
   double W = grid.get("width");	
   double xkl, ykl, zetaInitial, psiInitial;
   
   // set IC for zeta and psi (interior points)
   for (int l = 1 ; l <= ny - 2 ; l++)
   { 
    for (int k = 1 ; k <= nx - 2 ; k++)	
    {
     xkl     = grid.get("x",k,l);
     ykl     = grid.get("y",k,l);

     zetaInitial = xkl*xkl + ykl*ykl;
     grid.set("zeta",k,l,zetaInitial);
     psiInitial = 0;
     grid.set("psi",k,l,psiInitial);


    }

   }
   // set BC
   grid.set_psi_bc("diri",0.0,0.0,0.0,0.0);
   grid.set_zeta_bc();



grid.solve(ti,tf,reynolds,cfl,tol,omega,dumpPeriod,"sor",filename);
*/




