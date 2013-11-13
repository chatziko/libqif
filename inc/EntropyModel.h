#ifndef _QIF_EntropyModel_h_
#define _QIF_EntropyModel_h_
/*
This file belongs to the LIBQIF library.
A Quantitative Information Flow C++ Toolkit Library.
Copyright (C) 2013  Universidad Nacional de Río Cuarto(National University of Río Cuarto).
Author: Martinelli Fernán - fmartinelli89@gmail.com - Universidad Nacional de Río Cuarto (Argentina)
LIBQIF Version: 1.0
Date: 12th Nov 2013 
========================================================================
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

=========================================================================
*/
#include <fstream>
#include "Prob.h"
#include "Channel.h"
#include "Gain.h"
//#include <typeinfo>

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
		
		virtual DoubleType vulnerability(Prob& pi)=0;
		
		virtual DoubleType cond_vulnerability(Prob& pi)=0;
		
		virtual DoubleType leakage(Prob& pi)=0;
		
		virtual DoubleType entropy(Prob& pi)=0;
		
		virtual DoubleType cond_entropy(Prob& pi)=0;
		
		virtual DoubleType capacity()=0;
		
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

		virtual const char* class_name() { return "EntropyModel";}
		
	protected:
		int plotter_flag;  /*!< This integer will be used as a flag to determine with which plotter draw.
			by default will be -1. \n
			0 : SciLab \n
			1 : GNU-Plot \n
			2 : MatLab \n
			3 : Maple  */
	protected:
		Channel* C;	
	protected:
		Gain* g;	
};
#endif
