#include "GLeakage.h"
GLeakage::GLeakage(Channel c,Gain gain)
{
	C=&c;
	g=&gain;
}	
GLeakage::~GLeakage()
{
	C->~Channel();
	g->~Gain();
}	
void * GLeakage::compare_over_prior(Channel other_channel)
{
	
}
	
void * GLeakage::compare_over_gain(Channel other_channel,Prob prior)
{
	
}

//-------------- declaring the theoric algoritmhs implementation

double GLeakage::vulnerability(Prob pi)
{
	int w;
	int x;
	double sum_x;
	double max_w=0;
	for(w=0;w<g->guesses_number();w++)
	{
		sum_x=0;
		for(x=0;x<g->inputs_number();x++)
		{
			sum_x+= pi.at(x) * g->at(w,x);
		}
		if(sum_x>max_w){
			max_w=sum_x;
		}
	}
	return max_w;
}
		
double GLeakage::cond_vulnerability(Prob pi)
{
	int w;
	int x;
	int y;
	double sum_x;
	double max_w;
	double sum_y=0;
		
	for(y=0;y<C->outputs_number();y++)
	{
		max_w=0;
		for(w=0;w<g->guesses_number();w++)
		{
			sum_x=0;
			for(x=0;x<g->inputs_number();x++)
			{
				sum_x+= pi.at(x) * g->at(w,x) * C->at(x,y);
			}
			if(sum_x>max_w){
				max_w=sum_x;
			}
		}
		sum_y+= max_w;
	}
	return sum_y;
}	
	
double GLeakage::leakage(Prob pi)
{
	return log (cond_vulnerability(pi)/vulnerability(pi));
}	

double GLeakage::additive_leakage(Prob pi)
{
	return (entropy(pi)-cond_entropy(pi));
}
	
double GLeakage::entropy(Prob pi)
{
	return -log(vulnerability(pi));
}		
double GLeakage::cond_entropy(Prob pi)
{
	return -log(cond_vulnerability(pi));
}		
double GLeakage::capacity()
{
	throw 1; //It is not supported
}