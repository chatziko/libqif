#ifndef _QIF_EntropyModel_h_
#define _QIF_EntropyModel_h_

#include "Prob.h"

/*! \class EntropyModel
 *  \brief A generic model of entropy that defines the basic function to compute and plot.
 *
 *  For most information about the SciLab plotter engine see <a href="../papers/p1.pdf">here</a> \n
 *  For most information about the GNU-Plot plotter engine see <a href="../papers/p1.pdf">here</a> \n
 *  For most information about the MatLab plotter engine see <a href="../papers/p1.pdf">here</a> \n
 *  For most information about the Maple plotter engine see <a href="../papers/p1.pdf">here</a> \n
 */
class EntropyModel
{
	public:	
		//! A normal constructor. You will need choose a Engine before plotting. The SciLab plotter engine is choosen by default.
		/*!
		\sa ~EntropyModel(),change_to_GNUPlot(),change_to_MatLab(),change_to_Maple().
		*/
		EntropyModel();
		
		//! A normal destroyer member. This function closes the current engine.
		/*!
		\sa EntropyModel()
		*/
		~EntropyModel();
		
		//---------------------------------------------------------
		//theorethic algorithms
		
		virtual double vulnerability(Prob pi) = 0;
		
		virtual double cond_vulnerability(Prob pi) = 0;
		
		virtual double leakage(Prob pi)= 0;
		
		virtual double entropy(Prob pi)= 0;
		
		virtual double cond_entropy(Prob pi)= 0;
		
		virtual double capacity()= 0;
		
		//----------------------------------------------------------
		//plotter functions
		
		//! Plots the vulnerability function in 2D using the chosen engine.
		void plot2d_vulnerability();
		
		//! Plots the conditional vulnerability function in 2D using the chosen engine.
		void plot2d_cond_vulnerability();
		
		//! Plots the leakage function in 2D using the chosen engine.
		void plot2d_leakage();
		
		//! Plots the entropy function in 2D using the chosen engine.
		void plot2d_entropy();
		
		//! Plots the conditional entropy function in 2D using the chosen engine.
		void plot2d_cond_entropy();
		
		//! Plots the vulnerability function in 3D using the chosen engine.
		void plot3d_vulnerability();
		
		//! Plots the conditional vulnerability function in 3D using the chosen engine.
		void plot3d_cond_vulnerability();
		
		//! Plots the leakage function in 3D using the chosen engine.
		void plot3d_leakage();
		
		//! Plots the entropy function in 3D using the chosen engine.
		void plot3d_entropy();
		
		//! Plots the conditional entropy function in 3D using the chosen engine.
		void plot3d_cond_entropy();
		
		//--------------------------------------------------------
		
		//! Closes the previous engine, and opens the SciLab Engine.
		/*!
		\sa change_to_gnuplot(),change_to_matlab(),change_to_maple().
		*/
		void change_to_scilab();
		
		//! Closes the previous engine, and opens the GNU-Plot Engine.
		/*!
		\sa change_to_scilab(),change_to_matlab(),change_to_maple().
		*/
		void change_to_gnuplot();
		
		//! Closes the previous engine, and opens the MatLab Engine.
		/*!
		\sa change_to_scilab(),change_to_gnuplot(),change_to_maple().
		*/
		void change_to_matlab();
		
		//! Closes the previous engine, and opens the Maple Engine.
		/*!
		\sa change_to_scilab(),change_to_gnuplot(),change_to_matlab().
		*/
		void change_to_maple();
		
	protected:
		int plotter;  /*!< This integer will be used as a flag to determine with which plotter draw.
			by default will be -1. \n
			0 : SciLab \n
			1 : GNU-Plot \n
			2 : MatLab \n
			3 : Maple  */
};

#endif
