#include <iostream>
#include <fstream>
#include <cmath> 
#include <assert.h> 
#include <string>
#include "FluidElement.h"
using namespace std;

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

void FluidElement::set_rhs(double rhsSpec)
{
	rhs = rhsSpec;
}

void FluidElement::set_psiManu(double psiManuSpec)
{
	psiManu = psiManuSpec;
}


// CLEAR FUNCTION
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FluidElement::clear()
{

	psi  = 0;
	zeta = 0;
	rhs  = 0;
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

double FluidElement::get_rhs()
{
	return rhs;
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

