#include <iostream>
#include <fstream>
#include <cmath> 
#include <assert.h> 
#include <string>

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
Grid::Grid(int nxSpec,int nySpec, double lSpec, double wSpec)
{

	// ensure specified grid doesn't exceed allocated memory 
	assert(max(nx,ny) <= 200);
	
	dims = 2;
	nx = nxSpec;
	ny = nySpec;
	nz = 1;

	iterations = 0;

	l  = lSpec;
	w  = wSpec;

	dx = l/(nx - 1);
	dy = w/(ny - 1);
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


else if(var.compare("rhs") == 0)	
{
	nodes[k][l].rhs = val;	
}

else if(var.compare("psiManu") == 0)	
{
	nodes[k][l].psiManu = val;	
}

else 
{
	cout << "Invalid get_var selection (select: psi,zeta,rhs, psiManu) "  << 	endl;
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


else if(var.compare("rhs") == 0)	
{
	val = nodes[k][l].rhs;	
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
	cout << "Invalid get_var selection (select: psi,zeta,rhs,psiManu,x,y) "  << 	endl;
	val = -9999.999;
}

return val;

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

else 
{
	cout << "Invalid get_var selection (select: length,width,dx,dy) "  << 	endl;
	val = -9999.999;
}

return val;

}


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

errorPsi[0] = errorPsiMax;
errorPsi[1] = errorPsiCum/n;
}





// SOLVE FOR PSI IN POISSON EQUATION 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Grid::solve_poisson(double tol,double omega)
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
	double rhs  = 0;


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

			psikl   = Grid::get("psi",k,l);	
			psikp1l = Grid::get("psi",k+1,l);
			psikm1l = Grid::get("psi",k-1,l);	
			psiklp1 = Grid::get("psi",k,l+1);
			psiklm1 = Grid::get("psi",k,l-1);
			rhs  = Grid::get("rhs",k,l);

			//psi_new = dx*dx*dy*dy/(2*(dx*dx + dy*dy))*((psikp1l + psikm1l)/(dx*dx) + (psiklp1 + psiklm1)/(dy*dy) - rhs);
			psi_new = ( (psikp1l + psikm1l) + (psiklp1 + psiklm1)*beta*beta - rhs*dx*dx)/(1 + beta*beta)/2;
			residualPsiI = abs(psikl - psi_new); // residualPsi should acually be calculated after all psi_new calc
				
			if (residualPsiI > residualPsiMax)
			{
				residualPsiMax = residualPsiI;
			}
				
			//cout << residualPsiI << endl;
			psi_new = (1-omega)*psikl + omega*psi_new;
			Grid::set("psi",k,l,psi_new);
			}		
	
		}
	}
	residualPsi[0] = residualPsiMax;
	residualPsi[1] = residualPsiCum/iterations;
	Grid::calc_psi_error();
	
}


// CREATE HEADER FOR GRID NEW GRID PARAMETER DATA FILE
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::create_grid_params_header(std::string filename)
{
	ofstream myfile;	
	myfile.open (filename.c_str(),ios::app);
	myfile << "Length " << "Width " << "NX " << "NY " << "DX "  << "DY " << "ResidMax " << "ResidAvg " << "ErrorPsiMax " << "ErrorPsiAvg " << "DX*DX " << "Iterations " << endl;
	myfile.close();
}


// SAVE DATA TO GRID PARAMETER DATA FILE 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Grid::save_grid_params(std::string filename)
{
	ofstream myfile;
	myfile.open (filename.c_str(),ios::app);
	myfile << " " << l << " " << w << " " << nx << " " << ny << " " << dx << " " << dy << " " << residualPsi[0] << " " << residualPsi[1] << " " << errorPsi[0] << " " << errorPsi[1] << " " << dx*dx << " " << iterations << endl;
 	myfile.close();

/*
	ofstream myfile;
	
	if (myfile.is_open())
	{	

	}

	else
	{

	ofstream myfile;
	myfile.open (filename.c_str(),ios::app);	
	myfile << "Length " << "Width " << "NX " << "NY " << "DX "  << "DY " << "DX*DX " << "ResidMax " << "ResidAvg " << " " << endl;
	myfile << " " << l << " " << w << " " << nx << " " << ny << " " << dx << " " << dy << " " << dx*dx << " " << residualPsi[0] << " " << residualPsi[0] << " " << endl;
	myfile.close();
	}
*/
}


// Save vtk structured grid and node data
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void Grid::save_vtk_scaler_field(std::string filename)
{

assert (nx == ny); // edit this function to allow rectalinear grids.

    vtkSmartPointer<vtkImageData> imageData =
            vtkSmartPointer<vtkImageData>::New();

    imageData->SetDimensions(nx, ny, nz);

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

  int i = 0;
  for (int k = 0 ; k <= nx - 1 ; k++)
  {
	for (int l = 0 ; l <= ny - 1 ; l++)
	{
	double xI    = nodes[k][l].pos[0];
	double yI    = nodes[k][l].pos[1];
	double zI    = nodes[k][l].pos[2];
        double psiI = nodes[k][l].psi;
   	double psiManuI = nodes[k][l].psiManu;
	double zetaI = nodes[k][l].psi - nodes[k][l].psiManu;

        director->SetTuple3(i, xI, yI, zI);
        psiArray->SetValue(i, psiI);
	psiManuArray->SetValue(i, psiManuI);
	zetaArray->SetValue(i, zetaI);
	i ++;	
	}
  }


    imageData->GetPointData()->AddArray(director);
    director->SetName("Director");

    imageData->GetPointData()->AddArray(psiArray);
    psiArray->SetName("PSI");

    imageData->GetPointData()->AddArray(psiManuArray);
    psiManuArray->SetName("PSI_MANU");

    imageData->GetPointData()->AddArray(zetaArray);
    zetaArray->SetName("ZETA");



    vtkSmartPointer<vtkXMLImageDataWriter> writer =
            vtkSmartPointer<vtkXMLImageDataWriter>::New();

    writer->SetFileName(filename.c_str());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(imageData->GetProducerPort());
#else
    writer->SetInputData(imageData);
#endif
    writer->Write();
}
//vtk.1045678.n5.nabble.com/Proper-Usage-of-the-vtkStructuredGrid-td1241910.html
// the first part commented out is an attempt to use dx not equal dy
/*
// Declare structured grid and cooresponding points
  vtkSmartPointer<vtkStructuredGrid> structuredGrid =  
   vtkSmartPointer<vtkStructuredGrid>::New();
  
  vtkSmartPointer<vtkPoints> points = 
   vtkSmartPointer<vtkPoints>::New();

// Declare scalar feilds 
  vtkSmartPointer<vtkDoubleArray> scalar =  
   vtkSmartPointer<vtkDoubleArray>::New();
    scalar->SetNumberOfComponents(1);
    scalar->SetNumberOfTuples(nx * ny * nz);


  int counter = 0;
  for (int k = 0 ; k <= nx - 1 ; k++)
  {
	for (int l = 0 ; l <= ny - 1 ; l++)
	{
	double x = nodes[k][l].pos[0];
	double y = nodes[k][l].pos[1];
	double z = nodes[k][l].pos[2];
	points->InsertNextPoint(x, y, z);
	scalar->SetValue(counter, 1.45345);	
	counter ++;	
	}
  }


  structuredGrid->SetDimensions(nx,ny,nz);
  structuredGrid->SetPoints(points);

  structuredGrid.GetPointData().AddArray(scalar);
	
  vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  writer->SetFileName("test3.vtk");
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(structuredGrid);
#else
  writer->SetInputData(structuredGrid);
#endif
  writer->Write();

}




// Create a grid
  vtkSmartPointer<vtkStructuredGrid> structuredGrid =
    vtkSmartPointer<vtkStructuredGrid>::New();

  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkDoubleArray> energy = vtkSmartPointer<vtkDoubleArray>::New();

    energy->SetNumberOfComponents(1);
    energy->SetNumberOfTuples(nx * ny * nz);

//   scalars = vtk.vtkFloatArray()
//   vtkSmartPointer<vtkDoubleArray> energy =
//   vtkSmartPointer<vtkDoubleArray>::New();

int counter = 0;
for (int k = 0 ; k <= nx - 1 ; k++)
{
	for (int l = 0 ; l <= ny - 1 ; l++)
	{
	counter ++;
	double x = nodes[k][l].pos[0];
	double y = nodes[k][l].pos[1];
	double z = nodes[k][l].pos[2];
	points->InsertNextPoint(x, y, z);
	double e = counter;
	energy->SetValue(counter, e);	
	}
}


  structuredGrid->SetDimensions(nx,ny,nz);
  structuredGrid->SetPoints(points);

//    structuredGrid->GetPointData()->AddArray(energy);
//structuredGrid->GetPointData().SetScalars(energy);
//energy->SetName("PSI");

vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();

  //Write file

 // vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
 //   vtkSmartPointer<vtkXMLStructuredGridWriter>::New();

  writer->SetFileName("test2.vtk");
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(structuredGrid);
#else
  writer->SetInputData(structuredGrid);
#endif
  writer->Write();

}

*/

/*

{
    int nx = 10, ny = 10, nz = 10;

    vtkSmartPointer<vtkImageData> imageData =
            vtkSmartPointer<vtkImageData>::New();

    imageData->SetDimensions(nx, ny, nz);



    vtkSmartPointer<vtkDoubleArray> director =
            vtkSmartPointer<vtkDoubleArray>::New();

    director->SetNumberOfComponents(3);
    director->SetNumberOfTuples(nx * ny * nz);

    vtkSmartPointer<vtkDoubleArray> energy =
            vtkSmartPointer<vtkDoubleArray>::New();

    energy->SetNumberOfComponents(1);
    energy->SetNumberOfTuples(nx * ny * nz);

    for (int i = 0; i < director->GetNumberOfTuples(); ++i) {
        double t = 1.0;
        double p = 0.0;
        double e = 5.0;
        double x = sin(t) * cos(p),
                y = sin(t) * sin(p),
                z = cos(t);

        director->SetTuple3(i, x, y, z);
        energy->SetValue(i, e);
    }

    imageData->GetPointData()->AddArray(director);
    director->SetName("Director");

    imageData->GetPointData()->AddArray(energy);
    energy->SetName("Energy");

    vtkSmartPointer<vtkXMLImageDataWriter> writer =
            vtkSmartPointer<vtkXMLImageDataWriter>::New();

    writer->SetFileName("test.vti");
#if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(imageData->GetProducerPort());
#else
    writer->SetInputData(imageData);
#endif
    writer->Write();

//    return EXIT_SUCCESS;

*/

/*
vtkSmartPointer<vtkDoubleArray> director = vtkSmartPointer<vtkDoubleArray>::New();
director->SetNumberOfComponents(dims);
director->SetNumberOfTuples(nx*ny*nz);


vtkSmartPointer<vtkDoubleArray> psiVtk = vtkSmartPointer<vtkDoubleArray>::New();
psiVtk->SetNumberOfComponents(1);
psiVtk->SetNumberOfTuples(nx*ny*nz);

vtkSmartPointer<vtkDoubleArray> zetaVtk = vtkSmartPointer<vtkDoubleArray>::New();
zetaVtk->SetNumberOfComponents(1);
zetaVtk->SetNumberOfTuples(nx*ny*nz);




int counter = 0;
//    for (int i = 0; i < director->GetNumberOfTuples(); ++i) {
for (int k = 0 ; k <= nx - 1 ; k++)
{
	for (int l = 0 ; l <= ny - 1 ; l++)
	{

	counter ++;
        double x = nodes[k][l].pos[0];
        double y = nodes[k][l].pos[1];
        double psiI = nodes[k][l].psi;
 	double zetaI = nodes[k][l].zeta;
        director->SetTuple3(counter, x, y, 1);
        psiVtk->SetValue(counter, psiI);
	zetaVtk->SetValue(counter, zetaI);    
	}

}


director->GetPointData()->AddArray(psiVtk);
psiVtk->SetName("PSI");

director->GetPointData()->AddArray(zetaVtk);
zetaVtk->SetName("ZETA");








vtkSmartPointer<vtkXMLStructuredDataWriter> writer = vtkSmartPointer<vtkXMLStructuredDataWriter>::New();

    writer->SetFileName("test.vti");
#if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(StructuredData->GetProducerPort());
#else
    writer->SetInputData(StructuredData);
#endif
    writer->Write();

    return EXIT_SUCCESS;

*/
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Grid::set_psi_bc(double psiLbc, double psiRbc, double psiTbc, double psiBbc)
{
// set left and right
for (int l = 0 ; l <= ny - 1 ; l++)
 {

 nodes[0][l].psi    = psiLbc;
 nodes[nx-1][l].psi = psiRbc;
 }	

// set top and bottom
 for (int k = 0 ; k <= nx - 1 ; k++)
 {
 nodes[k][0].psi    = psiBbc;
 nodes[k][ny-1].psi = psiTbc;
 }
}





