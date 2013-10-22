#include "Shannon.h"

Shannon::Shannon(Channel c)
{
	C=&c;	
}
	
Shannon::~Shannon()
{
	C->~Channel();	
}
//-------------- declaring the theoric algoritmhs implementation


double Shannon::vulnerability(Prob pi)
{
	throw 1; //It is not supported
}
		
double Shannon::cond_vulnerability(Prob pi)
{
	throw 1; //It is not supported
}			

double Shannon::entropy(Prob pi)
{
	int x;
	double sum_x=0;
	for(x=0;x<C->inputs_number();x++){
		sum_x+= pi.at(x) * (log (pi.at(x)) / log(2));
	}
	return -sum_x;
	// - sum x pi(x) log2 p(x)
	//log2 p(x) = log p(x) / log (2)
}		

double Shannon::cond_entropy(Prob pi)
{
	int y,x;
	double sum_y=0;
	double sum_x;
	
	for(y=0;y<C->outputs_number();y++){
		sum_x=0;
		for(x=0; x<C->inputs_number();x++){
			sum_x+= C->at(x,y) * (log (C->at(x,y)) / log(2));
		}
		sum_y+= pi.at(y) - sum_x;
	}
	return sum_y;
	//sum y p(y) - sum x C[x,y] log2 C[x,y]
	//log2 C[x,y] = log C[x,y] / log (2)
}		

double Shannon::leakage(Prob pi)
{
	return (entropy(pi)-cond_entropy(pi));
}	

double Shannon::capacity()
{
	//
	//
	return 0;
}
