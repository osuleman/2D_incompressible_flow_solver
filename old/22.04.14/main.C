/*      
	FCFD Semester Project  
	Part 1  
	Omar Suleman 
	13/03/13  
*/

#include <iostream>
#include <fstream>
#include <cmath> 
using namespace std;



int main()
{

	ofstream myfile;
	myfile.open ("param_study.txt");

	
	myfile << " iteration_tolerance " << " SOR_param " << " nx " << " ny " << " dx*dx " << " iterations " << " average error " << endl;

	myfile.close();


for (int tol = 6 ; tol <= 6 ; tol++)
{

	for (double SOR = 9 ; SOR <= 9 ; SOR++)
	{

		for (int grid = 1 ; grid <= 10 ; grid++)
		{





// Universal Constants	

	double pi = atan(1)*4;
	cout << pi << endl;

// Define Geometry of Domain
	double L  = 1;
	double x_0 = 0;

	double W  = 1;
	double y_0 = 0;

	
// Define Grid Parameters 
	int nx = 10*grid;
	int ny = 10*grid;
	
	double dx   = L/(nx - 1);
	double dy   = W/(ny - 1);	
	double beta = dx/dy;
	
// Fix grid size and only use part of grid.
	double x             [200][200] = {};
	double y             [200][200] = {};
	double phi           [200][200] = {};
	double phi_fake      [200][200] = {};
	double RHS           [200][200] = {};
	
	
// Define BC
	double bottom_bc = 0;
	double top_bc    = 0;		
	double left_bc   = 0;
	double right_bc  = 0;


// Iteration Parameters
	double error_tol = pow(10,-tol);
	double omega = SOR/10;



// Initial Conditions (phi), made up solution (phi_fake), and analytical RHS of phi_fake (RHS) 	
	for (int k = 0 ; k <= nx - 1 ; k++)
	{
		for (int l = 0 ; l <= ny - 1 ; l++)
		{
			x[k][l]   = x_0 + k*dx;
			y[k][l]   = y_0 + l*dy;
			phi[k][l] = cos(x[k][l]*2*pi/L)*cos(y[k][l]*2*pi/W);
			phi_fake[k][l] = sin(x[k][l]*2*pi/L)*sin(y[k][l]*2*pi/W); // this is a fake solution to Phi
			RHS[k][l] = -4*pi*pi*sin(x[k][l]*2*pi/L)*sin(y[k][l]*2*pi/W)*(1/L/L + 1/W/W);
		}
	}





// Boundary Conditions 
	for (int l = 0 ; l <= nx - 1 ; l++)
	{
		phi[0][l] = left_bc;
		phi[nx-1][l] = right_bc;
	}

	for (int k = 0 ; k <= nx - 1 ; k++)
	{
		phi[k][0] = bottom_bc;
		phi[k][ny-1] = top_bc;
	}


/*
// Calculate RHS for phi_fake  (Laplacian(phi_fake))
	for (int k = 1 ; k <= nx - 2 ; k++)
	{
		for (int l = 1 ; l <= ny - 2 ; l++)
		{
		
		RHS[k][l] = (phi_fake[k+1][l] - 2*phi_fake[k][l] + phi_fake[k-1][l])/(dx*dx) + (phi_fake[k][l+1] - 2*phi_fake[k][l] + phi_fake[k][l-1])/(dy*dy);
					
		}
	}	
*/	

// Initialize Iteration Variables
	double phi_new 	[200][200] = {};
	double residual          = error_tol;
	double residual_max 	 = error_tol;




// make an iteration counter to kill loop if counter > max iterations
	
	int iterations = 0;
	while (residual_max >= error_tol)	
	{	
		// reset residual to find new max		
		residual_max = 0;
		iterations += 1;
		// cout << "iterations " << iterations << endl;
		for (int k = 1 ; k <= nx - 2 ; k++)
		{
			for (int l = 1 ; l <= ny - 2 ; l++)
			{
		
				phi_new[k][l] = ( (phi[k+1][l] + phi[k-1][l]) + (phi[k][l + 1] + phi[k][l - 1])*beta*beta - RHS[k][l]*dx*dx)/(1 + beta*beta)/2;
				
				// change residual to be cumulative error r = r + max(r(i)(j))
				
				residual = abs(phi[k][l] - phi_new[k][l]); // residual should acually be calculated after all phi_new calc
				
				if (residual > residual_max)
				{
					residual_max = residual;
				}
				
				phi[k][l] = (1-omega)*phi[k][l] + omega*phi_new[k][l];
				
				//cout << residual << endl;
			}
		}	
	}



// Calculate the absolute average error of solution
		double error_total = 0;
		double counter =0;
		for (int k = 1 ; k <= nx - 2 ; k++)
		{
			for (int l = 1 ; l <= ny - 2 ; l++)
			{


				error_total += abs(phi[k][l] - phi_fake[k][l]);
				//cout << "error i " << abs(RHS_Numerical[k][l] - RHS[k][l]) <<  "sum of error " << error_total << endl;
				counter += 1;
				
			}
		}

	cout << " counter" << counter << endl;
	double error_avg = error_total/counter;
	cout << "average error " << error_avg << endl;
	cout << "iterations " << iterations << endl;

// output results


	ofstream myfile;
	myfile.open ("param_study.txt",ios::app);

	
	myfile << " " << error_tol << " " << omega << " " << nx << " " << ny << " " << dx*dx << " " << iterations << " " << error_avg << " " << endl;

	myfile.close();


	ofstream writephi;
	writephi.open ("phi.txt",ios::app);

	writephi << endl << endl << endl << " iteration_tolerance " << error_tol <<  " SOR_param " << omega << " nx " << nx << " ny " << ny << " iterations " << iterations << " average error " << error_avg;	

	
	
		for (int k = 1 ; k <= ny - 2 ; k++)
		{
			writephi << endl;
			
			for (int l = 1 ; l <= nx - 2 ; l++)
			{
				writephi << phi[k][l] << " ";

				
			}
		}



	writephi.close();


		}
	}
}



return 0;
}





