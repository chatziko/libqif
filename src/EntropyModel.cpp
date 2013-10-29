#include "EntropyModel.h"
//#include <stack-c.h>
//#include <call_scilab.h>
/* Author: Fernan Martinelli*/
EntropyModel::EntropyModel()
{
	//Its not implemented yet
	plotter= -1;
}
		
EntropyModel::~EntropyModel()
{
	switch (plotter)
	{
		case -1: //do nothing
		break;
		
		case 0: 
		//close the scilab engine
//		if ( TerminateScilab(NULL) == FALSE ) {
//			throw 1;
//		}
		break;
		
		case 1: 
		//Its not implemented yet
		throw 1;
		break;
		
		case 2: 
		//Its not implemented yet
		throw 1;
		break;		
		
		case 3: 
		//Its not implemented yet
		throw 1;
		break;
	}
}
void EntropyModel::plot2d_vulnerability()
{
	switch (plotter)
	{
		case -1:
		//Must choose an engine before use it
		throw 1; 
		break;
		
		case 0: 
		//Its not implemented yet****************************
		break;
		
		case 1: 
		//Its not implemented yet
		throw 1;
		break;
		
		case 2: 
		//Its not implemented yet
		throw 1;
		break;		
		
		case 3: 
		//Its not implemented yet
		throw 1;
		break;
	}
}
void EntropyModel::plot2d_cond_vulnerability()
{
	switch (plotter)
	{
		case -1:
		//Must choose an engine before use it
		throw 1; 
		break;
		
		case 0: 
		//Its not implemented yet****************************
		break;
		
		case 1: 
		//Its not implemented yet
		throw 1;
		break;
		
		case 2: 
		//Its not implemented yet
		throw 1;
		break;		
		
		case 3: 
		//Its not implemented yet
		throw 1;
		break;
	}
}
void EntropyModel::plot2d_leakage()
{
	switch (plotter)
	{
		case -1:
		//Must choose an engine before use it
		throw 1; 
		break;
		
		case 0: 
		//Its not implemented yet****************************
		break;
		
		case 1: 
		//Its not implemented yet
		throw 1;
		break;
		
		case 2: 
		//Its not implemented yet
		throw 1;
		break;		
		
		case 3: 
		//Its not implemented yet
		throw 1;
		break;
	}
}
void EntropyModel::plot2d_entropy()
{
	switch (plotter)
	{
		case -1:
		//Must choose an engine before use it
		throw 1; 
		break;
		
		case 0: 
		//Its not implemented yet****************************
		break;
		
		case 1: 
		//Its not implemented yet
		throw 1;
		break;
		
		case 2: 
		//Its not implemented yet
		throw 1;
		break;		
		
		case 3: 
		//Its not implemented yet
		throw 1;
		break;
	}
}
void EntropyModel::plot2d_cond_entropy()
{
	switch (plotter)
	{
		case -1:
		//Must choose an engine before use it
		throw 1; 
		break;
		
		case 0: 
		//Its not implemented yet****************************
		break;
		
		case 1: 
		//Its not implemented yet
		throw 1;
		break;
		
		case 2: 
		//Its not implemented yet
		throw 1;
		break;		
		
		case 3: 
		//Its not implemented yet
		throw 1;
		break;
	}
}
void EntropyModel::plot3d_vulnerability()
{
	switch (plotter)
	{
		case -1:
		//Must choose an engine before use it
		throw 1; 
		break;
		
		case 0: 
		//Its not implemented yet****************************
		break;
		
		case 1: 
		//Its not implemented yet
		throw 1;
		break;
		
		case 2: 
		//Its not implemented yet
		throw 1;
		break;		
		
		case 3: 
		//Its not implemented yet
		throw 1;
		break;
	}
}
void EntropyModel::plot3d_cond_vulnerability()
{
	switch (plotter)
	{
		case -1:
		//Must choose an engine before use it
		throw 1; 
		break;
		
		case 0: 
		//Its not implemented yet****************************
		break;
		
		case 1: 
		//Its not implemented yet
		throw 1;
		break;
		
		case 2: 
		//Its not implemented yet
		throw 1;
		break;		
		
		case 3: 
		//Its not implemented yet
		throw 1;
		break;
	}
}
void EntropyModel::plot3d_leakage()
{
	switch (plotter)
	{
		case -1:
		//Must choose an engine before use it
		throw 1; 
		break;
		
		case 0: 
		//Its not implemented yet****************************
		break;
		
		case 1: 
		//Its not implemented yet
		throw 1;
		break;
		
		case 2: 
		//Its not implemented yet
		throw 1;
		break;		
		
		case 3: 
		//Its not implemented yet
		throw 1;
		break;
	}
}
void EntropyModel::plot3d_entropy()
{
	switch (plotter)
	{
		case -1:
		//Must choose an engine before use it
		throw 1; 
		break;
		
		case 0: 
		//Its not implemented yet****************************
		break;
		
		case 1: 
		//Its not implemented yet
		throw 1;
		break;
		
		case 2: 
		//Its not implemented yet
		throw 1;
		break;		
		
		case 3: 
		//Its not implemented yet
		throw 1;
		break;
	}
}
void EntropyModel::plot3d_cond_entropy()
{
	switch (plotter)
	{
		case -1:
		//Must choose an engine before use it
		throw 1; 
		break;
		
		case 0: 
		//Its not implemented yet****************************
		break;
		
		case 1: 
		//Its not implemented yet
		throw 1;
		break;
		
		case 2: 
		//Its not implemented yet
		throw 1;
		break;		
		
		case 3: 
		//Its not implemented yet
		throw 1;
		break;
	}
}
void EntropyModel::change_to_scilab()
{
	//Its not implemented yet
	switch (plotter)
	{
		case -1: plotter=0;
		//open scilab
//			#ifdef _MSC_VER
//			if ( StartScilab ( NULL, NULL, NULL ) == FALSE )
//			#else
//			if ( StartScilab ( getenv ( "SCI" ), NULL, NULL ) == FALSE )
//			#endif
//			{
//				throw 1;
//			}
		break;
		
		case 0: //do nothing
		break;
		
		case 1: //close gnuplot
		plotter= -1;
		change_to_scilab();
		break;
		
		case 2: //close matlab
		plotter= -1;
		change_to_scilab();
		break;		
		
		case 3: //close maple
		plotter= -1;
		change_to_scilab();
		break;
	}
}
void EntropyModel::change_to_gnuplot()
{
	//Its not implemented yet
	throw 1;
}		
void EntropyModel::change_to_matlab()
{
	//Its not implemented yet
	throw 1;
}
void EntropyModel::change_to_maple()
{
	//Its not implemented yet
	throw 1;
}