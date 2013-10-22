#ifndef _QIF_Mechanism_h_
#define _QIF_Mechanism_h_

#include "Channel.h"
#include "Graph.h"

class Mechanism : public Channel
{
	public:
		Mechanism(std::string& new_channel_elements,Graph new_graph);
		
		~Mechanism();
		
		bool is_diffenrential_private(double epsilon);
				
	protected:
		Graph graph;
};
#endif