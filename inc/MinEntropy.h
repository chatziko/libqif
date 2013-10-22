#ifndef _QIF_MinEntropy_h_
#define _QIF_MinEntropy_h_

#include "Prob.h"
#include "Channel.h"
#include "EntropyModel.h"
/*! \class MinEntropy
 *  \brief The Min-Entropy model of entropy.
 *
 *  For most information about the foundations of this theory see <a href="../papers/p1.pdf">here</a> 
 */
class MinEntropy : public EntropyModel
{
	public:
		
		MinEntropy(Channel c);
		
		~MinEntropy();
		
		double vulnerability(Prob pi);
			
		double cond_vulnerability(Prob pi);
			
		double leakage(Prob pi);
			
		double entropy(Prob pi);
			
		double cond_entropy(Prob pi);
			
		double capacity();
		
	protected:
		Channel* C;				
};
#endif
