#ifndef _QIF_Guessing_h_
#define _QIF_Guessing_h_

#include "Prob.h"
#include "Channel.h"
#include "EntropyModel.h"
/*! \class Guessing
 *  \brief The guessing model of entropy.
 *
 *  For most information about this theory see 
 */
class Guessing : public EntropyModel
{
	public:
		Guessing(Channel c);
		
		~Guessing();
		
		double vulnerability(Prob pi);
			
		double cond_vulnerability(Prob pi);
			
		double leakage(Prob pi);
			
		double entropy(Prob pi);
			
		double cond_entropy(Prob pi);
			
		double capacity();
		
	protected:
		//Channel C;		
};

#endif
