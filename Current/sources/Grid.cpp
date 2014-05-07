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

using namespace std;


// Constructor for 2D grid
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Grid::Grid()
{
	// set iterations to zero, and length to 1
	iterations = 0;
	l  = 1;
	w  = 1;
}



// Constructor for 2D grid
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Grid::Grid(int nxSpec,int nySpec, double lSpec, double wSpec, double lRefSpec, double uRefSpec)
{

	// ensure specified grid doesn't exceed allocated memory 
	assert(max(nx,ny) <= nMax);
	
	dims = 2;
	nx = nxSpec;
	ny = nySpec;
	nz = 1;

	iterations = 0;

	l  = lSpec;
	w  = wSpec;
	lRef = lRefSpec;
	uRef = uRefSpec;

	dx = l/(nx - 1);
	dy = w/(ny - 1);
	dz = 0;

	

	for (int k = 0 ; k <= nx - 1 ; k++)
	{
		for (int l = 0 ; l <= ny - 1 ; l++)
		{
			nodes[k][l].clear();
			double x = k*dx/lRef;
			double y = l*dy/lRef;
			double z = 0;
			nodes[k][l].set_position(x,y,z);
		}
	}

}

// Constructor for 2D grid
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::set_grid(int nxSpec,int nySpec, double lSpec, double wSpec, double lRefSpec, double uRefSpec)
{

	// ensure specified grid doesn't exceed allocated memory 
	assert(max(nx,ny) <= nMax);
	
	dims = 2;
	nx = nxSpec;
	ny = nySpec;
	nz = 1;

	iterations = 0;

	l  = lSpec;
	w  = wSpec;
	lRef = lRefSpec;
	uRef = uRefSpec;


// non dimensional dx and dy
	dx = l/lRef/(nx - 1);
	dy = w/lRef/(ny - 1);
	dz = 0;

	

	for (int k = 0 ; k <= nx - 1 ; k++)
	{
		for (int l = 0 ; l <= ny - 1 ; l++)
		{
			nodes[k][l].clear();
			double x = k*dx;
			double y = l*dy;
			double z = 0;
			nodes[k][l].set_position(x, y,z);
		}
	}

}
// GET FluidElement POINTER
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
FluidElement *Grid::get_node(int k, int l)
{
	return &nodes[k][l];
}


// SET FluidElement VARIABLE AT SPECIFIED POS IN NODES TO SPECIFIED VAL
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::set(std::string var, int k, int l, double val)
{

if(var.compare("psi") == 0)	
{
	nodes[k][l].psi = val;
}


else if(var.compare("zeta") == 0)	
{
	nodes[k][l].zeta = val;
}


else if(var.compare("rhsPsi") == 0)	
{
	nodes[k][l].rhsPsi = val;	
}

else if(var.compare("psiManu") == 0)	
{
	nodes[k][l].psiManu = val;	
}

else 
{
	cout << "Invalid get_var selection (select: psi,zeta,rhsPsi, psiManu) "  << 	endl;
	exit(EXIT_FAILURE);
}

}



// GET FluidElement VARIABLE AT SPECIFIED POSITION IN NODES 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid::get(std::string var, int k, int l)
{


double val;

if(var.compare("psi") == 0)	
{
	val = nodes[k][l].psi;
}


else if(var.compare("zeta") == 0)	
{
	val = nodes[k][l].zeta;
}


else if(var.compare("rhsPsi") == 0)	
{
	val = nodes[k][l].rhsPsi;	
}

else if(var.compare("psiManu") == 0)	
{
	val = nodes[k][l].psiManu;	
}

else if(var.compare("x") == 0)	
{
	val = nodes[k][l].pos[0];	
}

else if(var.compare("y") == 0)	
{
	val = nodes[k][l].pos[1];	
}

else 
{
	cout << "Invalid get_var selection (select: psi,zeta,rhsPsi,psiManu,x,y) "  << 	endl;
	val = -9999.999;
	exit(EXIT_FAILURE);
}

return val;

}


// set Grid PARAMETER 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::set(std::string var, double val)
 {

 if(var.compare("length") == 0)		
 {
	l = val;	
 }

 else if(var.compare("width") == 0)	
 {
	w = val;	
 }

 else if(var.compare("t") == 0)	
 {
	val = t;	
 }
 else 
 {
	cout << "Invalid get_var selection (select: length,width,t,uLeft,uRight,uTop,Ubottom) "  << 	endl;
	exit(EXIT_FAILURE);
 }


}



// GET Grid PARAMETER 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid::get(std::string var)
{


double val;

if(var.compare("length") == 0)		
{
	val = l;	
}

else if(var.compare("width") == 0)	
{
	val = w;	
}

else if(var.compare("dx") == 0)	
{
	val = dx;	
}

else if(var.compare("dy") == 0)	
{
	val = dy;	
}

else if(var.compare("iterations") == 0)	
{
	val = iterations;	
}

else if(var.compare("t") == 0)	
{
	val = t;	
}

else if(var.compare("errorPsiMax") == 0)	
{
	val = errorPsi[0];	
}

else if(var.compare("errorPsiAvg") == 0)	
{
	val = errorPsi[1];
}


else 
{
	cout << "Invalid get_var selection (select: length,width,dx,dy, iterations, t) "  << 	endl;
	val = -9999.999;
	exit(EXIT_FAILURE);
}

return val;

}

//SET BC FOR PSI (only diri BC set-up)  Takes dimensional bc and non dimensionalizes them.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::set_psi_bc(std::string type,double psiLbc, double psiRbc, double psiTbc, double psiBbc)
{

 if(type.compare("diri") == 0)	
 {
  for (int l = 0 ; l <= ny - 1 ; l++)
  {
   // set left and right
   nodes[0][l].psi    = psiLbc/uRef/lRef;
   nodes[nx-1][l].psi = psiRbc/uRef/lRef;
  }	

  // set top and bottom
  for (int k = 0 ; k <= nx - 1 ; k++)
  {
   nodes[k][0].psi    = psiBbc/uRef/lRef;
   nodes[k][ny-1].psi = psiTbc/uRef/lRef;
  }	
 }

 else 
 {
  cout << "select diri boundary condition type" << endl;
  exit(EXIT_FAILURE);
 }
}



//SET BC FOR ZETA (only diri BC set-up)  Takes dimensional bc and non dimensionalizes them.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::set_zeta_bc(std::string type, double zetaLbc, double zetaRbc, double zetaTbc, double zetaBbc)
{

 if(type.compare("diri") == 0)	
 {
  for (int l = 0 ; l <= ny - 1 ; l++)
  {
   // set left and right
   nodes[0][l].zeta    = zetaLbc*lRef/uRef;
   nodes[nx-1][l].zeta = zetaRbc*lRef/uRef;
  }	

  // set top and bottom
  for (int k = 0 ; k <= nx - 1 ; k++)
  {
   nodes[k][0].zeta    = zetaBbc*lRef/uRef;
   nodes[k][ny-1].zeta = zetaTbc*lRef/uRef;
  }	
 }


 else 
 {
  cout << "select diri or vel boundary condition type" << endl;
  exit(EXIT_FAILURE);
 }
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::set_boundary_velocities(double vLeft, double vRight, double uTop, double uBottom)
{
 // set top and bottom
 for (int k = 1 ; k <= nx - 2 ; k++)
 {
  // set top boundary velocity
  nodes[k][0].u    = uBottom/uRef;
  nodes[k][0].v    = 0;

  // set bottom boundary velocity
  nodes[k][ny-1].u = uTop/uRef;
  nodes[k][ny-1].v = 0;
 }
 for (int l = 0 ; l <= ny - 1 ; l++)
 {
  // set left boundary velocity
  nodes[0][l].u = 0;
  nodes[0][l].v = vLeft/uRef;

  // set right boundary velocity
  nodes[nx-1][l].u = 0;
  nodes[nx-1][l].v = vRight/uRef;
 }	
}



//CALC PSI ERROR
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::calc_psi_error()
{
// initialize errorPsi and counter
double errorPsiCum = 0.0;
double errorPsiI = 0.0;
double errorPsiMax = 0.0;
int n = 0;

	for (int k = 1 ; k <= nx - 2 ; k++)
	{
		for (int l = 1 ; l <= ny - 2 ; l++)
		{
		n += 1;
		errorPsiI = abs(nodes[k][l].psi - nodes[k][l].psiManu);

		//cout << nodes[k][l].pos[0] << " " << nodes[k][l].pos[1] << " " << errorPsiI << endl; 
		if (errorPsiI > errorPsiMax)
		 	{
			errorPsiMax = errorPsiI;
			}
		errorPsiCum += errorPsiI;
		}
	}


// store grid params with dimensions (dimensionalize before returning)
errorPsi[0] = errorPsiMax*uRef*lRef;
errorPsi[1] = errorPsiCum/n*uRef*lRef;
}

//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::set_psi_bc()
{  
 // set top and bottom
 for (int k = 0 ; k <= nx - 1 ; k++)
 {
  nodes[k][0].psi    = nodes[k][0].u*dx + nodes[k][1].psi;
  nodes[k][ny-1].psi = nodes[k][ny-1].u*dx + nodes[k][ny-2].psi;
 }
 for (int l = 1 ; l <= ny - 2 ; l++)
 {
  // set left and right
  nodes[0][l].psi    = nodes[0][l].v*dy + nodes[1][l].psi;
  nodes[nx-1][l].psi = nodes[nx-1][l].v*dy + nodes[nx-2][l].psi;
 }	
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::set_zeta_bc()
{
  // set top and bottom
  for (int k = 0 ; k <= nx - 1 ; k++)
  {
   nodes[k][0].zeta    =  2*(nodes[k][0].psi - nodes[k][1].psi + dy*nodes[k][0].u/uRef)/dy/dy;
   nodes[k][ny-1].zeta =  2*(nodes[k][ny-1].psi - nodes[k][ny-2].psi - dy*nodes[k][ny-1].u/uRef)/dy/dy;
  }


  for (int l = 1; l <= ny - 2; l++)
  {
   // set left and right
   nodes[0][l].zeta    =  2*(nodes[0][l].psi - nodes[1][l].psi - dx*nodes[0][l].v/uRef)/dx/dx;
   nodes[nx-1][l].zeta =  2*(nodes[nx-1][l].psi - nodes[nx-2][l].psi + dx*nodes[nx-1][l].v/uRef)/dx/dx;
  }	


}





// TEST POISSON EQUATION (solve poisson equation with rhs(psi_manu)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid::test_poisson(double tol,double omega)
{

// initialize residualPsi and iterations
	double residualPsiMax = tol;
	double residualPsiCum = 0;
	double residualPsiI   = 0;
	iterations = 0;

	

// Initialize Stensil Variables
	double beta = dx/dy;	
	double psi_new = 0;
	double psikl   = 0;	
	double psikp1l = 0;
	double psikm1l = 0;	
	double psiklp1 = 0;
	double psiklm1 = 0;
	double rhsPsi  = 0;


	while (residualPsiMax >= tol)	
	{	
		// reset residualPsi to find new max
		residualPsiMax = 0;		
		iterations += 1;
		
		// loop through interior points
		for (int k = 1 ; k <= nx - 2 ; k++)
		{
			for (int l = 1 ; l <= ny - 2 ; l++)
			{

			psikl   = nodes[k][l].psi;	
			psikp1l = nodes[k+1][l].psi;
			psikm1l = nodes[k-1][l].psi;	
			psiklp1 = nodes[k][l+1].psi;
			psiklm1 = nodes[k][l-1].psi;
			rhsPsi     = nodes[k][l].rhsPsi;

			psi_new = ( (psikp1l + psikm1l) + (psiklp1 + psiklm1)*beta*beta - rhsPsi*dx*dx)/(1 + beta*beta)/2;
			residualPsiI = abs(psikl - psi_new); // residualPsi should acually be calculated after all psi_new calc
				
			if (residualPsiI > residualPsiMax)
			{
				residualPsiMax = residualPsiI;
			}
				
			//cout << residualPsiI << endl;
			psi_new = (1-omega)*psikl + omega*psi_new;
			nodes[k][l].psi = psi_new;
			}		
	
		}
	}
	residualPsi[0] = residualPsiMax*uRef*lRef;
	residualPsi[1] = residualPsiCum/iterations*uRef*lRef;
	Grid::calc_psi_error();
	
}





// SOLVE FOR PSI IN POISSON EQUATION 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid::solve_poisson(double tol,double omega)
{
http://www.nhl.com/
// initialize residualPsi and iterations
	double residualPsiMax = tol;
	double residualPsiCum = 0;
	double residualPsiI   = 0;
	iterations = 0;

	

// Initialize Stensil Variables
	double beta = dx/dy;	
	double psi_new = 0;
	double psikl   = 0;	
	double psikp1l = 0;
	double psikm1l = 0;	
	double psiklp1 = 0;
	double psiklm1 = 0;
	double rhs     = 0;


	while (residualPsiMax >= tol)	
	{	
		// reset residualPsi to find new max
		residualPsiMax = 0;		
		iterations += 1;
		
		// loop through interior points
		for (int k = 1 ; k <= nx - 2 ; k++)
		{
			for (int l = 1 ; l <= ny - 2 ; l++)
			{

			psikl   = nodes[k][l].psi;	
			psikp1l = nodes[k+1][l].psi;
			psikm1l = nodes[k-1][l].psi;	
			psiklp1 = nodes[k][l+1].psi;
			psiklm1 = nodes[k][l-1].psi;
			rhs     =-nodes[k][l].zeta;

			psi_new = ( (psikp1l + psikm1l) + (psiklp1 + psiklm1)*beta*beta - rhs*dx*dx)/(1 + beta*beta)/2;
			residualPsiI = abs(psikl - psi_new); // residualPsi should acually be calculated after all psi_new calc
				
			if (residualPsiI > residualPsiMax)
			{
				residualPsiMax = residualPsiI;
			}
				
			//cout << residualPsiI << endl;
			psi_new = (1-omega)*psikl + omega*psi_new;
			nodes[k][l].psi = psi_new;
			}		
	
		}
	}
	residualPsi[0] = residualPsiMax*uRef*lRef;
	residualPsi[1] = residualPsiCum/iterations*uRef*lRef;
	Grid::calc_psi_error();
	
}



// calculate 2nd order accurate, centered 1st or 2nd spatial derivative 
// Only works if DX and DY Constant
double Grid::diff(int k, int l, std::string f, std::string  xi, int deriv) 
{
 // first derivative
 if (deriv == 1)
 {
	// x direction
	if(xi.compare("x") == 0)		
	{
	double kp1l  = Grid::get(f,k+1,l);
	double km1l  = Grid::get(f,k-1,l);
	double ans = (kp1l - km1l)/2/dx;
	}

	// y direction
	else if(xi.compare("y") == 0)	
	{
	double klp1  = Grid::get(f,k,l+1);
	double klm1  = Grid::get(f,k,l-1);
	double ans = (klp1 - klm1)/2/dy;		
	}
	else
	{
	cout << "invalid spatial derivative.  choose x or y direction " << endl;
	exit(EXIT_FAILURE);
	}



 }
 else if (deriv == 2)
 {

	// x direction
	if(xi.compare("x") == 0)		
	{
	double kl    = Grid::get(f,k,l);	
	double kp1l  = Grid::get(f,k+1,l);
	double km1l  = Grid::get(f,k-1,l);
	double ans = (kp1l - 2*kl + km1l)/dx/dx;
	}

	// y direction
	else if(xi.compare("y") == 0)	
	{
	double kl    = Grid::get(f,k,l);	
	double klp1  = Grid::get(f,k,l+1);
	double klm1  = Grid::get(f,k,l-1);
	double ans = (klp1 - 2*kl + klm1)/dy/dy;		
	}
	else
	{
	cout << "invalid spatial derivative.  choose x or y direction " << endl;
	exit(EXIT_FAILURE);
	}	



}
else 
{

cout << "invalid derivative selection.  Choose 1st or 2nd order only" << endl;
exit(EXIT_FAILURE);

}


}




// SOLVE NON DIM. VORTICITY EQUATION FOR DZETA_DT
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Grid::solve_dzeta_dt(int k, int l, double reynolds)
{

nodes[k][l].u =  Grid::diff(k,l,"psi","y",1);
nodes[k][l].v = -Grid::diff(k,l,"psi","x",1);
//cout << "u" << nodes[k][l].u << " v" << nodes[k][l].v << endl;



// diffusive terms
double ans = (Grid::diff(k,l,"zeta","x",2) + Grid::diff(k,l,"zeta","y",2))/reynolds;

// convective terms
ans = ans - nodes[k][l].u*Grid::diff(k,l,"zeta","x",1) - nodes[k][l].v*Grid::diff(k,l,"zeta","y",1); 


return ans;
}



// Runge-Kutta 4 Solver
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::rk4_vorticity(double dt, double reynolds)
{

FluidElement nodesOld[nMax][nMax];
//nodesOld = nodes;
//Grid gridOld = grid;
/*
double zetaOld = nodes[k][l].zeta;
// double tOld = t;
double rhs = Grid::solve_dzeta_dt(int k, int l, double reynolds);
nodes[k][l].zeta = zetaOld + 0.5*dt*rhs;
double rhsSum = rhs;

// t = told + dt/2

rhs = Grid::solve_dzeta_dt(int k, int l, double reynolds);
nodes[k][l].zeta = zetaOld + 0.5*dt*rhs;
rhsSum = rhsSum + 2.0*rhs;
rhs = Grid::solve_dzeta_dt(int k, int l, double reynolds);
nodes[k][l].zeta = zetaOld + dt*rhs;
rhsSum = rhsSum + 2.0*rhs;

// t = told + dt 

rhs = Grid::solve_dzeta_dt(int k, int l, double reynolds);
rhsSum = rhsSum + rhs;
nodes[k][l].zeta = zetaOld + rhsSum/6.0;

}


double zetaT = Grid::get("zeta",k,l);
double ans = zetaT + Grid::solve_dzeta_dt(k,l,visc);
t = t + dt;
return ans;
*/
}


// SOLVE NON DIM. VORTICITY EQUATION FOR ZETA(t+1)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::solve_zeta_tp1(double dt, double reynolds)
{
 for (int k = 1 ; k <= nx - 2; k++)
 {
	for (int l = 1 ; l <= ny - 2 ; l++)
	{
	 double zetaT = nodes[k][l].zeta;
	 double ans = zetaT + Grid::solve_dzeta_dt(k,l,reynolds)*dt;
	 nodes[k][l].zetaTp1 = ans;
	}
 }


// save zetaTp1 to zeta
 for (int k = 1 ; k <= nx - 2; k++)
 {
	for (int l = 1 ; l <= ny - 2 ; l++)
	{
	 nodes[k][l].zeta = nodes[k][l].zetaTp1;
	}
 }


// increment time
 t = t + dt;


}







// SOLVE VORTICITY EQUATION FOR ZETA(t+1)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::solve(double ti, double tf, double reynolds, double cfl, double tol, double omega, double dumpPeriod, std::string filenameSpec)
{
//Grid::solve_poisson(tol,omega);
Grid::set_zeta_bc();
t = ti;
double counter = 0;
int save    = 0;
double dt   = 0;
 if (dx >= dy)
 {
  dt = cfl*dx/uRef;
 }
else
 {
  dt = cfl*dy/uRef;
 }

// create filename for t=0 and save initial grid
  std::ostringstream os;
  os << filenameSpec << counter << ".vti";
  // "_dt:" << dt << 
  std::string filename = os.str();
  Grid::save_vtk_field(filename);


while (t <= tf)
{
 counter ++;
 Grid::solve_zeta_tp1(dt,reynolds);
 Grid::solve_poisson(tol,omega);
 Grid::set_zeta_bc();
 save = fmod(counter,dumpPeriod); 
 if(save == 0)
 {
  std::ostringstream os;
  os << filenameSpec << counter << ".vti";
  filename = os.str();
  Grid::save_vtk_field(filename);
//  Grid::save_txt_field(filename);
 }
 
}



}




// CREATE HEADER FOR GRID NEW GRID PARAMETER DATA FILE
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::create_grid_params_header(std::string filenameSpec)
{
	// generate filename
	std::ostringstream os;
	os << filenameSpec << "_L:" << l << "_W:" << l << ".txt";
	std::string filename = os.str();

	ofstream myfile;	
	myfile.open (filename.c_str(),ios::app);
	myfile << "NX " << "NY " << "DX "  << "DY " << "ResidMax " << "ResidAvg " << "ErrorPsiMax " << "ErrorPsiAvg " << "DX*DX " << "Iterations " << endl;
	myfile.close();
}


// SAVE DATA TO GRID PARAMETER DATA FILE non dim data
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::save_grid_params(std::string filenameSpec)
{
	// generate filename
	std::ostringstream os;
	os << filenameSpec << "_L:" << l << "_W:" << w << ".txt";
	std::string filename = os.str();


	ofstream myfile;
	myfile.open (filename.c_str(),ios::app);
	myfile << nx << " " << ny << " " << dx << " " << dy << " " << residualPsi[0] << " " << residualPsi[1] << " " << errorPsi[0] << " " << errorPsi[1] << " " << dx*dx << " " << iterations << endl;
 	myfile.close();

}





void Grid::save_txt_field(std::string filenameSpec)
{
	// generate filename
	std::ostringstream os;
	os << filenameSpec << "_vectorData"<< "_nx:" << nx << "_ny:" << ny << ".txt";
	std::string filename = os.str();

	ofstream myfile;
	myfile.open (filename.c_str(),ios::app);

	for (int k = 1 ; k <= ny - 2 ; k++)
	{
		myfile << endl;
		
		for (int l = 1 ; l <= nx - 2 ; l++)
		{
			myfile << nodes[k][l].u << " ";
				
		}
	}
}


// Save vtk structured grid and node data non dim data
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Grid::save_vtk_field(std::string filenameSpec)
{

assert (nx == ny); // edit this function to allow rectalinear grids.

    vtkSmartPointer<vtkImageData> scalarData =
            vtkSmartPointer<vtkImageData>::New();

    scalarData->SetDimensions(nx, ny, nz);

    vtkSmartPointer<vtkImageData> vectorData =
            vtkSmartPointer<vtkImageData>::New();

    vectorData->SetDimensions(nx, ny, nz);

    vtkSmartPointer<vtkDoubleArray> director =
            vtkSmartPointer<vtkDoubleArray>::New();

    director->SetNumberOfComponents(3);
    director->SetNumberOfTuples(nx * ny * nz);


    vtkSmartPointer<vtkDoubleArray> psiArray =
            vtkSmartPointer<vtkDoubleArray>::New();

    psiArray->SetNumberOfComponents(1);
    psiArray->SetNumberOfTuples(nx * ny * nz);


    vtkSmartPointer<vtkDoubleArray> psiManuArray =
            vtkSmartPointer<vtkDoubleArray>::New();

    psiManuArray->SetNumberOfComponents(1);
    psiManuArray->SetNumberOfTuples(nx * ny * nz);


    vtkSmartPointer<vtkDoubleArray> zetaArray =
            vtkSmartPointer<vtkDoubleArray>::New();

    zetaArray->SetNumberOfComponents(1);
    zetaArray->SetNumberOfTuples(nx * ny * nz);


    vtkSmartPointer<vtkDoubleArray> velocityArray =
            vtkSmartPointer<vtkDoubleArray>::New();

    velocityArray->SetNumberOfComponents(3);
    velocityArray->SetNumberOfTuples(nx * ny * nz);


  int i = 0;
  for (int k = 0 ; k <= nx - 1 ; k++)
  {
	for (int l = 0 ; l <= ny - 1 ; l++)
	{
	double xI       = nodes[k][l].pos[0];
	double yI       = nodes[k][l].pos[1];
	double zI       = nodes[k][l].pos[2];
        double psiI     = nodes[k][l].psi;
   	double psiManuI = nodes[k][l].psiManu;
	double zetaI    = nodes[k][l].zeta;
	double uI       = nodes[k][l].v;
	double vI       = nodes[k][l].u;
	double wI       = nodes[k][l].w;
	

        director->SetTuple3(i, xI, yI, zI);
        psiArray->SetValue(i, psiI);
	psiManuArray->SetValue(i, psiManuI);
	zetaArray->SetValue(i, zetaI);
	velocityArray->SetTuple3(i,uI,vI,wI);	
	i ++;	
	}
  }


    scalarData->GetPointData()->AddArray(director);
    vectorData->GetPointData()->AddArray(director);
    director->SetName("POSITION");

    vectorData->GetPointData()->AddArray(velocityArray);
    velocityArray->SetName("VELOCITY");

    scalarData->GetPointData()->AddArray(psiArray);
    psiArray->SetName("PSI");

    scalarData->GetPointData()->AddArray(psiManuArray);
    psiManuArray->SetName("PSI_MANU");

    scalarData->GetPointData()->AddArray(zetaArray);
    zetaArray->SetName("ZETA");



	// generate filename
	std::ostringstream os;
	os << filenameSpec << "_scalarData"<< "_nx:" << nx << "_ny:" << ny << ".vti";
	std::string filename = os.str();




    vtkSmartPointer<vtkXMLImageDataWriter> writer =
            vtkSmartPointer<vtkXMLImageDataWriter>::New();

    writer->SetFileName(filename.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(scalarData->GetProducerPort());
#else
    writer->SetInputData(scalarData);
#endif
    writer->Write();

	// generate filename
	std::ostringstream os2;
	os2 << filenameSpec << "_vectorData"<< "_nx:" << nx << "_ny:" << ny << ".vti";
	filename = os2.str();




    vtkSmartPointer<vtkXMLImageDataWriter> writer2 =
            vtkSmartPointer<vtkXMLImageDataWriter>::New();

    writer2->SetFileName(filename.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer2->SetInputConnection(vectorData->GetProducerPort());
#else
    writer2->SetInputData(vectorData);
#endif
    writer2->Write();


}







