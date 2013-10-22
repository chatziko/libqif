#ifndef _QIF_Shannon_h_
#define _QIF_Shannon_h_

#include "Prob.h"
#include "Channel.h"
#include "EntropyModel.h"
/*! \class Shannon
 *  \brief The shannon model of entropy.
 *
 *  For most information about the foundations of this theory see <a href="../papers/p1.pdf">here</a> 
 */
class Shannon : public EntropyModel
{
	public:
		Shannon(Channel c);
		
		~Shannon();
		
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
