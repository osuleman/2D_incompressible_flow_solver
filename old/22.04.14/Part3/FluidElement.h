#include <iostream>
#include <fstream>
#include <cmath> 
using namespace std;



class FluidElement
{
protected:
	double pos[2];
//	double psi, zeta, rhs;



public:
	
	void set_position(double x, double y);
	
/*
	void clear();
	void set_psi(double psiSpec);
	void set_zeta(double zetaSpec); 
	void set_rhs(double rhsSpec);
*/
}





void FluidElement::set_position(double x, double y);
{
	pos[0] = x;
	pos[1] = y;
}

