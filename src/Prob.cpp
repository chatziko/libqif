#include "Prob.h"

Prob::Prob (double* new_vector )
{
	//VER COMO IMPLEMENTAR ESTO.. NO ESTA DIRECTAMENTE EN ARMADILLO
}

Prob::Prob (char* new_vector )
{
	V=vec(new_vector);
}

Prob::~Prob()
{
	V.~vec();
}

int Prob::size()
{
	return V.size();
}

double Prob::at ( int index )
{
	return V.at(index); //OJO : VER SI ES NECESARIO CONTROLAR EL RANGO.
}