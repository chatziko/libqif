#ifndef _QIF_LinearProgram_h_
#define _QIF_LinearProgram_h_

#include <armadillo>
using namespace arma;
#include <stdlib.h>
#include <glpk.h>


class LinearProgram
{
	public:
	
	vec solve(char* equality,char * inequality,char * objective);
	
	vec solve(mat equality,mat inequality,vec objective);
	
	vec solve(char* equality,char * inequality,char * objective,char * rows_constraints);
	
	vec solve(mat equality,mat inequality,vec objective,mat rows_constraints);
};

#endif
