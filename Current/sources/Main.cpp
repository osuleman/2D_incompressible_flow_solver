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
// Grid Refinement Study
/*
int nMin = 10;
int nGrids = 20;
int dn  = 10;
GridRefinementStudy grs(nMin,nGrids,dn);
grs.test_derivative1();
grs.test_derivative2();

// testing different SOR Parameters, omega

int    nOmega = 3;
double omegaI = 0.6;
double dOmega = 0.3;
double omega = 0;


for (int i = 0; i <= nOmega-1; i ++)
{
 double tol = pow(10,-10);
        omega = omegaI + i*dOmega;



//grs.test_poisson(omega,tol);
}

*/
//grid params
int nx        = 30;
int ny        = 30;
double length = 1;
double width  = 1;
double lRef   = length;

double ti = 0;
double tf = 10;
double reynolds = 10;
double cfl = 0.2;
double tol = pow(10,-10);
double omega = 1;
double dumpPeriod = 1;
std::string filename = "results/test";

double vLbc, vRbc, uTbc, uBbc;

//lid velocity
vLbc = 0;
vRbc = 0;
uTbc = 1;
uBbc = 0;
double uRef   = uTbc;


Grid grid(nx,ny,length,width,lRef,uRef);
grid.define_wall_velocities(vLbc,vRbc,uTbc,uBbc);

   double L = grid.get("length");
   double W = grid.get("width");	
   double xkl, ykl, zetaInitial, psiInitial;
   
   // set IC for zeta and psi
   for (int l = 1 ; l <= ny - 2 ; l++)
   { 
    for (int k = 1 ; k <= nx - 2 ; k++)	
    {
     xkl     = grid.get("x",k,l);
     ykl     = grid.get("y",k,l);

     zetaInitial = 0;
     grid.set("zeta",k,l,zetaInitial);
     psiInitial = 0;
     grid.set("psi",k,l,zetaInitial);


    }

   }
   // set BC
   grid.set_psi_bc();
   grid.set_zeta_bc();



grid.solve(ti,tf,reynolds,cfl,tol,omega,dumpPeriod,filename);

}



