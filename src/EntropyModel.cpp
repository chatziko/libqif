#include "EntropyModel.h"
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
EntropyModel::EntropyModel()
{
	plotter_flag=-1;
	/* This integer will be used as a flag to determine with which plotter draw.
			by default will be -1. \n
			0 : SciLab \n
			1 : GNU-Plot \n
			2 : MatLab \n
			3 : Maple  */
}
		
EntropyModel::~EntropyModel()
{
	switch (plotter_flag)
	{
		case -1: //do nothing
		break;
		
		case 0: //do nothing
		break;
		
		case 1: //Its not implemented yet
		break;
		 
		case 2: //Its not implemented yet
		break;		
		
		case 3: //Its not implemented yet
		break;
	}
}

void EntropyModel::plot2d_vulnerability()
{
	//correct size control
	if(C->inputs_number()!=2)
	{
		throw 1; // X must be equal for both
	}

	//Should check for the specific classes if can implement it or not.
	//Shannon and Guessing have not vulnerability
	std::string className;
	className=this->class_name();
	if(className.compare("Shannon")==0){ throw 1;}
	if(className.compare("Guessing")==0){ throw 1;}
	switch (plotter_flag)
	{
		case -1:
			//Must choose an engine before use it
			throw 1; 
		break;
		
		case 0: 
		{
			std::ofstream outputFile;
			std::ifstream inputFile;
			std::string fileName;

			// creating the output file
			fileName = className+"_plot2d_vulnerability.sci";
			outputFile.open(fileName.c_str(),std::ofstream::out | std::ofstream::app);

			// importing the vulnerability function to the file
			inputFile.open("../sci_files/2dvulnerability.sci",std::ifstream::in);
			std::string current_line;
			while(getline(inputFile,current_line))
			{
				outputFile << current_line << std::endl;
			}
			inputFile.close();
			//writing channel C
			std::string channel_c="C = [";
			channel_c +=C->str;
			channel_c +="];";
			outputFile << channel_c << std::endl;

			//if it is GLeakage we should write Gain Function g
			if(className.compare("GLeakage")==0)
			{ 
				std::string gain_g="g = [";
				gain_g += g->str;
				gain_g +="];";
				outputFile << gain_g << std::endl;
			}

			// generating the plot
			std::string plotting_command="[a,b]=size(pi);";
			outputFile << plotting_command << std::endl;

			plotting_command="Res=zeros(b,a);";
			outputFile << plotting_command << std::endl;

			plotting_command="for i=1:b";
			outputFile << plotting_command << std::endl;
			
			//choosing the correct vulnerability function
			if(className.compare("GLeakage")==0)
			{
				plotting_command="Res(i)=V(pi(i));";
				outputFile << plotting_command << std::endl;
			}else{
				//it is MinEntropy
				plotting_command="Res(i)=Vm(pi(i));";
				outputFile << plotting_command << std::endl;
			}

   			plotting_command="end";				
			outputFile << plotting_command << std::endl;

			plotting_command="scf();";
			outputFile << plotting_command << std::endl;

			plotting_command="plot(pi,Res);";
			outputFile << plotting_command << std::endl;
			// end of generating the plot

			outputFile.close();
		}
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
	//correct size control
	if(C->inputs_number()!=2)
	{
		throw 1; // X must be equal for both
	}

	//Should check for the specific classes if can implement it or not.
	//Shannon and Guessing have not conditional vulnerability
	std::string className;
	className=this->class_name();
	if(className.compare("Shannon")==0){ throw 1;}
	if(className.compare("Guessing")==0){ throw 1;}
	switch (plotter_flag)
	{
		case -1:
			//Must choose an engine before use it
			throw 1; 
		break;
		
		case 0: 
		{
			std::ofstream outputFile;
			std::ifstream inputFile;
			std::string fileName;

			// creating the output file
			fileName = className+"_plot2d_cond_vulnerability.sci";
			outputFile.open(fileName.c_str(),std::ofstream::out | std::ofstream::app);

			// importing the cond_vulnerability function to the file
			inputFile.open("../sci_files/2dcond_vulnerability.sci",std::ifstream::in);
			std::string current_line;
			while(getline(inputFile,current_line))
			{
				outputFile << current_line << std::endl;
			}
			inputFile.close();
			//writing channel C
			std::string channel_c="C = [";
			channel_c +=C->str;
			channel_c +="];";
			outputFile << channel_c << std::endl;

			//if it is GLeakage we should write Gain Function g
			if(className.compare("GLeakage")==0)
			{ 
				std::string gain_g="g = [";
				gain_g += g->str;
				gain_g +="];";
				outputFile << gain_g << std::endl;
			}

			// generating the plot
			std::string plotting_command="[a,b]=size(pi);";
			outputFile << plotting_command << std::endl;

			plotting_command="Res=zeros(b,a);";
			outputFile << plotting_command << std::endl;

			plotting_command="for i=1:b";
			outputFile << plotting_command << std::endl;
			
			//choosing the correct cond_vulnerability function
			if(className.compare("GLeakage")==0)
			{
				plotting_command="Res(i)=Vp(pi(i));";
				outputFile << plotting_command << std::endl;
			}else{
				//it is MinEntropy
				plotting_command="Res(i)=Vpm(pi(i));";
				outputFile << plotting_command << std::endl;
			}

   			plotting_command="end";				
			outputFile << plotting_command << std::endl;

			plotting_command="scf();";
			outputFile << plotting_command << std::endl;

			plotting_command="plot(pi,Res);";
			outputFile << plotting_command << std::endl;
			// end of generating the plot

			outputFile.close();
		}
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
	//correct size control
	if(C->inputs_number()!=2)
	{
		throw 1; // X must be equal for both
	}
	std::string className;
	className=this->class_name();
	switch (plotter_flag)
	{
		case -1:
			//Must choose an engine before use it
			throw 1; 
		break;
		
		case 0: 
		{
			std::ofstream outputFile;
			std::ifstream inputFile;
			std::string fileName;

			// creating the output file
			fileName = className+"_plot2d_leakage.sci";
			outputFile.open(fileName.c_str(),std::ofstream::out | std::ofstream::app);

			// importing the leakage function to the file
			inputFile.open("../sci_files/2dleakage.sci",std::ifstream::in);
			std::string current_line;
			while(getline(inputFile,current_line))
			{
				outputFile << current_line << std::endl;
			}
			inputFile.close();
			//writing channel C
			std::string channel_c="C = [";
			channel_c +=C->str;
			channel_c +="];";
			outputFile << channel_c << std::endl;

			//if it is GLeakage we should write Gain Function g
			if(className.compare("GLeakage")==0)
			{ 
				std::string gain_g="g = [";
				gain_g += g->str;
				gain_g +="];";
				outputFile << gain_g << std::endl;
			}

			// generating the plot
			std::string plotting_command="[a,b]=size(pi);";
			outputFile << plotting_command << std::endl;

			plotting_command="Res=zeros(b,a);";
			outputFile << plotting_command << std::endl;

			plotting_command="for i=1:b";
			outputFile << plotting_command << std::endl;
			
			//choosing the correct leakage function
			if(className.compare("GLeakage")==0)
			{
				plotting_command="Res(i)=ML(pi(i));";
				outputFile << plotting_command << std::endl;
			}else{
				if(className.compare("MinEntropy")==0){
					//it is MinEntropy
					plotting_command="Res(i)=MLm(pi(i));";
					outputFile << plotting_command << std::endl;
				}else{
					if(className.compare("Guessing")==0){
						//it is Guessing
						plotting_command="Res(i)=MLg(pi(i));";
						outputFile << plotting_command << std::endl;
					}else{
						//it is Shannon
						plotting_command="Res(i)=MLs(pi(i));";
						outputFile << plotting_command << std::endl;
					}
				}
			} // end if

   			plotting_command="end";				
			outputFile << plotting_command << std::endl;

			plotting_command="scf();";
			outputFile << plotting_command << std::endl;

			plotting_command="plot(pi,Res);";
			outputFile << plotting_command << std::endl;
			// end of generating the plot

			outputFile.close();
		}
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
	//correct size control
	if(C->inputs_number()!=2)
	{
		throw 1; // X must be equal for both
	}
	std::string className;
	className=this->class_name();
	switch (plotter_flag)
	{
		case -1:
			//Must choose an engine before use it
			throw 1; 
		break;
		
		case 0: 
		{
			std::ofstream outputFile;
			std::ifstream inputFile;
			std::string fileName;

			// creating the output file
			fileName = className+"_plot2d_entropy.sci";
			outputFile.open(fileName.c_str(),std::ofstream::out | std::ofstream::app);

			// importing the entropy function to the file
			inputFile.open("../sci_files/2dentropy.sci",std::ifstream::in);
			std::string current_line;
			while(getline(inputFile,current_line))
			{
				outputFile << current_line << std::endl;
			}
			inputFile.close();
			//writing channel C
			std::string channel_c="C = [";
			channel_c +=C->str;
			channel_c +="];";
			outputFile << channel_c << std::endl;

			//if it is GLeakage we should write Gain Function g
			if(className.compare("GLeakage")==0)
			{ 
				std::string gain_g="g = [";
				gain_g += g->str;
				gain_g +="];";
				outputFile << gain_g << std::endl;
			}

			// generating the plot
			std::string plotting_command="[a,b]=size(pi);";
			outputFile << plotting_command << std::endl;

			plotting_command="Res=zeros(b,a);";
			outputFile << plotting_command << std::endl;

			plotting_command="for i=1:b";
			outputFile << plotting_command << std::endl;
			
			//choosing the correct entropy function
			if(className.compare("GLeakage")==0)
			{
				plotting_command="Res(i)=E(pi(i));";
				outputFile << plotting_command << std::endl;
			}else{
				if(className.compare("MinEntropy")==0){
					//it is MinEntropy
					plotting_command="Res(i)=Em(pi(i));";
					outputFile << plotting_command << std::endl;
				}else{
					if(className.compare("Guessing")==0){
						//it is Guessing
						plotting_command="Res(i)=Eg(pi(i));";
						outputFile << plotting_command << std::endl;
					}else{
						//it is Shannon
						plotting_command="Res(i)=Es(pi(i));";
						outputFile << plotting_command << std::endl;
					}
				}
			} // end if

   			plotting_command="end";				
			outputFile << plotting_command << std::endl;

			plotting_command="scf();";
			outputFile << plotting_command << std::endl;

			plotting_command="plot(pi,Res);";
			outputFile << plotting_command << std::endl;
			// end of generating the plot

			outputFile.close();
		}
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
	//correct size control
	if(C->inputs_number()!=2)
	{
		throw 1; // X must be equal for both
	}
	std::string className;
	className=this->class_name();
	switch (plotter_flag)
	{
		case -1:
			//Must choose an engine before use it
			throw 1; 
		break;
		
		case 0: 
		{
			std::ofstream outputFile;
			std::ifstream inputFile;
			std::string fileName;

			// creating the output file
			fileName = className+"_plot2d_cond_entropy.sci";
			outputFile.open(fileName.c_str(),std::ofstream::out | std::ofstream::app);

			// importing the cond_entropy function to the file
			inputFile.open("../sci_files/2dcond_entropy.sci",std::ifstream::in);
			std::string current_line;
			while(getline(inputFile,current_line))
			{
				outputFile << current_line << std::endl;
			}
			inputFile.close();
			//writing channel C
			std::string channel_c="C = [";
			channel_c +=C->str;
			channel_c +="];";
			outputFile << channel_c << std::endl;

			//if it is GLeakage we should write Gain Function g
			if(className.compare("GLeakage")==0)
			{ 
				std::string gain_g="g = [";
				gain_g += g->str;
				gain_g +="];";
				outputFile << gain_g << std::endl;
			}

			// generating the plot
			std::string plotting_command="[a,b]=size(pi);";
			outputFile << plotting_command << std::endl;

			plotting_command="Res=zeros(b,a);";
			outputFile << plotting_command << std::endl;

			plotting_command="for i=1:b";
			outputFile << plotting_command << std::endl;
			
			//choosing the correct cond_entropy function
			if(className.compare("GLeakage")==0)
			{
				plotting_command="Res(i)=Ep(pi(i));";
				outputFile << plotting_command << std::endl;
			}else{
				if(className.compare("MinEntropy")==0){
					//it is MinEntropy
					plotting_command="Res(i)=Epm(pi(i));";
					outputFile << plotting_command << std::endl;
				}else{
					if(className.compare("Guessing")==0){
						//it is Guessing
						plotting_command="Res(i)=Epg(pi(i));";
						outputFile << plotting_command << std::endl;
					}else{
						//it is Shannon
						plotting_command="Res(i)=Eps(pi(i));";
						outputFile << plotting_command << std::endl;
					}
				}
			} // end if

   			plotting_command="end";				
			outputFile << plotting_command << std::endl;

			plotting_command="scf();";
			outputFile << plotting_command << std::endl;

			plotting_command="plot(pi,Res);";
			outputFile << plotting_command << std::endl;
			// end of generating the plot

			outputFile.close();
		}
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
	//correct size control
	if(C->inputs_number()!=3)
	{
		throw 1; // X must be equal for both
	}

	std::string className;
	className=this->class_name();
	if(className.compare("Shannon")==0){ throw 1;}
	if(className.compare("Guessing")==0){ throw 1;}
	switch (plotter_flag)
	{
		case -1:
			//Must choose an engine before use it
			throw 1; 
		break;
		
		case 0: 
		{
			std::ofstream outputFile;
			std::ifstream inputFile;
			std::string fileName;

			// creating the output file
			fileName = className+"_plot3d_vulneability.sci";
			outputFile.open(fileName.c_str(),std::ofstream::out | std::ofstream::app);

			// importing the vulnerability function to the file
			inputFile.open("../sci_files/3dvulneability.sci",std::ifstream::in);
			std::string current_line;
			while(getline(inputFile,current_line))
			{
				outputFile << current_line << std::endl;
			}
			inputFile.close();
			//writing channel C
			std::string channel_c="C = [";
			channel_c +=C->str;
			channel_c +="];";
			outputFile << channel_c << std::endl;

			//if it is GLeakage we should write Gain Function g
			if(className.compare("GLeakage")==0)
			{ 
				std::string gain_g="g = [";
				gain_g += g->str;
				gain_g +="];";
				outputFile << gain_g << std::endl;
			}

			// generating the plot
			std::string plotting_command="scf();";
			outputFile << plotting_command << std::endl;

			//choosing the correct vulnerability function
			if(className.compare("GLeakage")==0)
			{
				plotting_command="fplot3d(pi,pi,V,30,45,'X@Y@Z');";
				outputFile << plotting_command << std::endl;
			}else{
				//it is MinEntropy
				plotting_command="fplot3d(pi,pi,Vm,30,45,'X@Y@Z');";
				outputFile << plotting_command << std::endl;
			} // end if			
			// end of generating the plot

			outputFile.close();
		}
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
	//correct size control
	if(C->inputs_number()!=3)
	{
		throw 1; // X must be equal for both
	}

	std::string className;
	className=this->class_name();
	if(className.compare("Shannon")==0){ throw 1;}
	if(className.compare("Guessing")==0){ throw 1;}
	switch (plotter_flag)
	{
		case -1:
			//Must choose an engine before use it
			throw 1; 
		break;
		
		case 0: 
		{
			std::ofstream outputFile;
			std::ifstream inputFile;
			std::string fileName;

			// creating the output file
			fileName = className+"_plot3d_cond_vulneability.sci";
			outputFile.open(fileName.c_str(),std::ofstream::out | std::ofstream::app);

			// importing the cond_vulnerability function to the file
			inputFile.open("../sci_files/3dcond_vulneability.sci",std::ifstream::in);
			std::string current_line;
			while(getline(inputFile,current_line))
			{
				outputFile << current_line << std::endl;
			}
			inputFile.close();
			//writing channel C
			std::string channel_c="C = [";
			channel_c +=C->str;
			channel_c +="];";
			outputFile << channel_c << std::endl;

			//if it is GLeakage we should write Gain Function g
			if(className.compare("GLeakage")==0)
			{ 
				std::string gain_g="g = [";
				gain_g += g->str;
				gain_g +="];";
				outputFile << gain_g << std::endl;
			}

			// generating the plot
			std::string plotting_command="scf();";
			outputFile << plotting_command << std::endl;

			//choosing the correct cond_vulnerability function
			if(className.compare("GLeakage")==0)
			{
				plotting_command="fplot3d(pi,pi,Vp,30,45,'X@Y@Z');";
				outputFile << plotting_command << std::endl;
			}else{
				//it is MinEntropy
				plotting_command="fplot3d(pi,pi,Vpm,30,45,'X@Y@Z');";
				outputFile << plotting_command << std::endl;
			} // end if			
			// end of generating the plot

			outputFile.close();
		}
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
	//correct size control
	if(C->inputs_number()!=3)
	{
		throw 1; // X must be equal for both
	}

	std::string className;
	className=this->class_name();
	switch (plotter_flag)
	{
		case -1:
			//Must choose an engine before use it
			throw 1; 
		break;
		
		case 0: 
		{
			std::ofstream outputFile;
			std::ifstream inputFile;
			std::string fileName;

			// creating the output file
			fileName = className+"_plot3d_leakage.sci";
			outputFile.open(fileName.c_str(),std::ofstream::out | std::ofstream::app);

			// importing the leakage function to the file
			inputFile.open("../sci_files/3dleakage.sci",std::ifstream::in);
			std::string current_line;
			while(getline(inputFile,current_line))
			{
				outputFile << current_line << std::endl;
			}
			inputFile.close();
			//writing channel C
			std::string channel_c="C = [";
			channel_c +=C->str;
			channel_c +="];";
			outputFile << channel_c << std::endl;

			//if it is GLeakage we should write Gain Function g
			if(className.compare("GLeakage")==0)
			{ 
				std::string gain_g="g = [";
				gain_g += g->str;
				gain_g +="];";
				outputFile << gain_g << std::endl;
			}

			// generating the plot
			std::string plotting_command="scf();";
			outputFile << plotting_command << std::endl;

			//choosing the correct leakage function
			if(className.compare("GLeakage")==0)
			{
				plotting_command="fplot3d(pi,pi,ML,30,45,'X@Y@Z');";
				outputFile << plotting_command << std::endl;
			}else{
				if(className.compare("MinEntropy")==0){
					//it is MinEntropy
					plotting_command="fplot3d(pi,pi,MLm,30,45,'X@Y@Z');";
					outputFile << plotting_command << std::endl;
				}else{
					if(className.compare("Guessing")==0){
						//it is Guessing
						plotting_command="fplot3d(pi,pi,MLg,30,45,'X@Y@Z');";
						outputFile << plotting_command << std::endl;
					}else{
						//it is Shannon
						plotting_command="fplot3d(pi,pi,MLs,30,45,'X@Y@Z');";
						outputFile << plotting_command << std::endl;
					}
				}
			} // end if			
			// end of generating the plot

			outputFile.close();
		}
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
	//correct size control
	if(C->inputs_number()!=3)
	{
		throw 1; // X must be equal for both
	}

	std::string className;
	className=this->class_name();
	switch (plotter_flag)
	{
		case -1:
			//Must choose an engine before use it
			throw 1; 
		break;
		
		case 0:
		{
			std::ofstream outputFile;
			std::ifstream inputFile;
			std::string fileName;

			// creating the output file
			fileName = className+"_plot3d_entropy.sci";
			outputFile.open(fileName.c_str(),std::ofstream::out | std::ofstream::app);

			// importing the entropy function to the file
			inputFile.open("../sci_files/3dentropy.sci",std::ifstream::in);
			std::string current_line;
			while(getline(inputFile,current_line))
			{
				outputFile << current_line << std::endl;
			}
			inputFile.close();
			//writing channel C
			std::string channel_c="C = [";
			channel_c +=C->str;
			channel_c +="];";
			outputFile << channel_c << std::endl;

			//if it is GLeakage we should write Gain Function g
			if(className.compare("GLeakage")==0)
			{ 
				std::string gain_g="g = [";
				gain_g += g->str;
				gain_g +="];";
				outputFile << gain_g << std::endl;
			}

			// generating the plot
			std::string plotting_command="scf();";
			outputFile << plotting_command << std::endl;

			//choosing the correct entropy function
			if(className.compare("GLeakage")==0)
			{
				plotting_command="fplot3d(pi,pi,E,30,45,'X@Y@Z');";
				outputFile << plotting_command << std::endl;
			}else{
				if(className.compare("MinEntropy")==0){
					//it is MinEntropy
					plotting_command="fplot3d(pi,pi,Em,30,45,'X@Y@Z');";
					outputFile << plotting_command << std::endl;
				}else{
					if(className.compare("Guessing")==0){
						//it is Guessing
						plotting_command="fplot3d(pi,pi,Eg,30,45,'X@Y@Z');";
						outputFile << plotting_command << std::endl;
					}else{
						//it is Shannon
						plotting_command="fplot3d(pi,pi,Es,30,45,'X@Y@Z');";
						outputFile << plotting_command << std::endl;
					}
				}
			} // end if			
			// end of generating the plot

			outputFile.close();
		}
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
	//correct size control
	if(C->inputs_number()!=3)
	{
		throw 1; // X must be equal for both
	}
	
	std::string className;
	className=this->class_name();
	switch (plotter_flag)
	{
		case -1:
			//Must choose an engine before use it
			throw 1; 
		break;
		
		case 0: 
		{
			std::ofstream outputFile;
			std::ifstream inputFile;
			std::string fileName;

			// creating the output file
			fileName = className+"_plot3d_cond_entropy.sci";
			outputFile.open(fileName.c_str(),std::ofstream::out | std::ofstream::app);

			// importing the cond_entropy function to the file
			inputFile.open("../sci_files/3dcond_entropy.sci",std::ifstream::in);
			std::string current_line;
			while(getline(inputFile,current_line))
			{
				outputFile << current_line << std::endl;
			}
			inputFile.close();
			//writing channel C
			std::string channel_c="C = [";
			channel_c +=C->str;
			channel_c +="];";
			outputFile << channel_c << std::endl;

			//if it is GLeakage we should write Gain Function g
			if(className.compare("GLeakage")==0)
			{ 
				std::string gain_g="g = [";
				gain_g += g->str;
				gain_g +="];";
				outputFile << gain_g << std::endl;
			}

			// generating the plot
			std::string plotting_command="scf();";
			outputFile << plotting_command << std::endl;

			//choosing the correct cond_entropy function
			if(className.compare("GLeakage")==0)
			{
				plotting_command="fplot3d(pi,pi,Ep,30,45,'X@Y@Z');";
				outputFile << plotting_command << std::endl;
			}else{
				if(className.compare("MinEntropy")==0){
					//it is MinEntropy
					plotting_command="fplot3d(pi,pi,Epm,30,45,'X@Y@Z');";
					outputFile << plotting_command << std::endl;
				}else{
					if(className.compare("Guessing")==0){
						//it is Guessing
						plotting_command="fplot3d(pi,pi,Epg,30,45,'X@Y@Z');";
						outputFile << plotting_command << std::endl;
					}else{
						//it is Shannon
						plotting_command="fplot3d(pi,pi,Eps,30,45,'X@Y@Z');";
						outputFile << plotting_command << std::endl;
					}
				}
			} // end if			
			// end of generating the plot

			outputFile.close();
		}
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
	switch (plotter_flag)
	{
		case -1: 
			plotter_flag=0;
			//open scilab
		break;
		
		case 0: //do nothing
		break;
		
		case 1: //close gnuplot
			plotter_flag=-1;
			change_to_scilab();
		break;
		
		case 2: //close matlab
			plotter_flag=-1;
			change_to_scilab();
		break;		
		
		case 3: //close maple
			plotter_flag=-1;
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