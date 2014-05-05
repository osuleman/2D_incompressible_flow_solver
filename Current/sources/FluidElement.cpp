#include <iostream>
#include <fstream>
#include <cmath> 
#include <assert.h> 
#include <string>
#include "FluidElement.h"
using namespace std;


FluidElement::FluidElement()
{

	psi      = 0;
	zeta     = 0;
	u        = 0;
	v        = 0;
	w        = 0;
	rhsPsi   = 0;
	psiManu  = 0;	
	pos[3]   = {};


}

/*
FluidElement* FluidElement::operator= (const FluidElement &equals)
{
return (*equals);

}
*/


// SET FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FluidElement::set_position(double x, double y, double z)
{
	pos[0] = x;
	pos[1] = y;
	pos[2] = z;
}

void FluidElement::set_psi(double psiSpec)
{
	psi = psiSpec;
}

void FluidElement::set_zeta(double zetaSpec)
{
	zeta = zetaSpec;
}

void FluidElement::set_rhsPsi(double rhsPsiSpec)
{
	rhsPsi = rhsPsiSpec;
}

void FluidElement::set_psiManu(double psiManuSpec)
{
	psiManu = psiManuSpec;
}

void FluidElement::set_u(double uSpec)
{
	u = uSpec;
}

void FluidElement::set_v(double vSpec)
{
	v = vSpec;
}

void FluidElement::set_w(double wSpec)
{
	w = wSpec;
}





// CLEAR FUNCTION
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FluidElement::clear()
{

	psi      = 0;
	zeta     = 0;
	u        = 0;
	v        = 0;
	w        = 0;
	rhsPsi   = 0;
	psiManu  = 0;
}


// GET FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double FluidElement::get_psi()
{
	return psi;
}

double FluidElement::get_zeta()
{
	return zeta;
}

double FluidElement::get_u()
{
	return u;
}

double FluidElement::get_v()
{
	return v;
}

double FluidElement::get_w()
{
	return w;
}

double FluidElement::get_rhsPsi()
{
	return rhsPsi;
}

double FluidElement::get_psiManu()
{
	return psiManu;
}

double FluidElement::get_x()
{
	return pos[0];
}

double FluidElement::get_y()
{
	return pos[1];
}

double FluidElement::get_z()
{
	return pos[2];
}

