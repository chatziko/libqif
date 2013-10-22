#include "MinEntropy.h"

MinEntropy::MinEntropy(Channel c)
{
	C=&c;
}
	
MinEntropy::~MinEntropy()
{
	C->~Channel();
}
//-------------- declaring the theoric algoritmhs implementation


double MinEntropy::vulnerability(Prob pi)
{
	int x;
	double max_x=0;
	for ( x = 0;x < C->inputs_number();x++ )
	{
		if ( pi.at(x) > max_x )
		{
			max_x = pi.at(x);
		}
	}
	return max_x;
}
		
double MinEntropy::cond_vulnerability(Prob pi)
{
	int x;
	int y;
	double sum_x;
	double max_x;
	double sum_y=0;
		
	for(y=0;y<C->outputs_number();y++)
	{
		max_x=0;
		for(x=0;x<C->inputs_number();x++)
		{
			sum_x= pi.at(x) * C->at(x,y);
			if(sum_x>max_x){
				max_x=sum_x;
			}
		}
		sum_y+= max_x;
	}
	return sum_y;
}		
double MinEntropy::leakage(Prob pi)
{
	return (entropy(pi)-cond_entropy(pi));
}		
double MinEntropy::entropy(Prob pi)
{
	return -log(vulnerability(pi)); 
}		
double MinEntropy::cond_entropy(Prob pi)
{
	return -log(cond_vulnerability(pi));
}		
double MinEntropy::capacity()
{
	int x;
	int y;
	double sum_x;
	double max_x;
	double sum_y=0;
		
	for(y=0;y<C->outputs_number();y++)
	{
		max_x=0;
		for(x=0;x<C->inputs_number();x++)
		{
			sum_x= C->at(x,y);
			if(sum_x>max_x){
				max_x=sum_x;
			}
		}
		sum_y+= max_x;
	}
	return log(sum_y);
}
