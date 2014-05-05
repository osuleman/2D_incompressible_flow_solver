#include <iostream>
#include <fstream>
#include <cmath> 
using namespace std;


class Grid
{
protected:

	const int nx, ny;
	const double dx, dy, l, w;
	FluidElement nodes[nx][ny];

public:

	Grid(int nxSpec, int nySpec)


}


Grid::Grid(int nxSpec, int nySpec, lSpec, wSpec)
{
	nx = nxSpec;
	ny = nySpec;

	l  = lSpec;
	w  = wSpec;

	dx   = L/(nx - 1);
	dy   = W/(ny - 1);

	for (int k = 0 ; k <= nx - 1 ; k++)
	{
		for (int l = 0 ; l <= ny - 1 ; l++)
		{
			
			double x = k*dx;
			double y = l*dy;
			nodes[k][l].set_position(x, y);


}

	 
	
