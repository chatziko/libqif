#ifndef _QIF_GLeakage_h_
#define _QIF_GLeakage_h_

#include "Prob.h"
#include "Channel.h"
#include "Gain.h"
#include "EntropyModel.h"
#include "LinearProgram.h"

/*! \class GLeakage
 *  \brief The Generalized Gain Function model of entropy.
 *
 *  For most information about the foundations of this theory see <a href="../papers/gleakage.pdf">here</a> 
 */
class GLeakage : public EntropyModel
{
	public:
		
		//! A normal constructor taking 2 arguments.
			/*!
			\sa ~GLeakage().
			*/
		GLeakage(Channel c,Gain g);
	
		//! A normal destroyer member.
			/*!
			\sa GLeakage()
			*/
		~GLeakage();
		
		double vulnerability(Prob pi);
			
		double cond_vulnerability(Prob pi);
			
		double leakage(Prob pi);
			
		double additive_leakage(Prob pi);
		
		double entropy(Prob pi);
			
		double cond_entropy(Prob pi);
			
		double capacity();
		
		void * compare_over_prior(Channel other_channel);
		
		void * compare_over_gain(Channel other_channel,Prob prior);
		
	protected:
		Gain* g;
		Channel* C;		
};

#endif
