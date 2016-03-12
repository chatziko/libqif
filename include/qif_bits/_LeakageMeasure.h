

// OBSOLETE
//
// Kostas:
// old base class, before restructuring the code in namespaces.
// Not used now, there might be useful code for plotting.



/*! \class LeakageMeasure
 *  \brief A generic model of entropy that defines the basic function to compute and plot.
 *
 *  For most information about the SciLab plotter engine see <a href="../papers/p1.pdf">here</a> \n
 *  For most information about the GNU-Plot plotter engine see <a href="../papers/p1.pdf">here</a> \n
 *  For most information about the MatLab plotter engine see <a href="../papers/p1.pdf">here</a> \n
 *  For most information about the Maple plotter engine see <a href="../papers/p1.pdf">here</a> \n
 */

template<typename eT>
class LeakageMeasure {
	public:
		Chan<eT> C;

		//! A normal constructor. You will need choose a Engine before plotting. The SciLab plotter engine is choosen by default.
		/*!
		\sa ~LeakageMeasure(),change_to_GNUPlot(),change_to_MatLab(),change_to_Maple().
		*/
		LeakageMeasure()					: C( ) {}
		LeakageMeasure(const Chan<eT>& C)	: C(C) {}
		LeakageMeasure(Chan<eT>&& C)		: C(C) {}

		//---------------------------------------------------------
		//theorethic algorithms

		virtual eT vulnerability(const Prob<eT>&) { throw std::runtime_error("not supported"); }

		virtual eT cond_vulnerability(const Prob<eT>&) { throw std::runtime_error("not supported"); }

		virtual eT entropy(const Prob<eT>& pi) = 0;

		virtual eT cond_entropy(const Prob<eT>& pi) = 0;

		virtual eT mult_leakage (const Prob<eT>& pi) { return cond_vulnerability(pi) / vulnerability(pi); }
		virtual eT leakage      (const Prob<eT>& pi) { return entropy(pi)            - cond_entropy(pi);  }

		virtual eT max_mult_leakage() { throw std::runtime_error("not supported"); }
		virtual eT capacity()         { throw std::runtime_error("not supported"); }

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

		virtual const char* class_name() {
			return "LeakageMeasure";
		}

		eT max_diff = def_max_diff<eT>();
		eT max_rel_diff = def_max_rel_diff<eT>();

	protected:
		int plotter_flag = -1;  /*!< This integer will be used as a flag to determine with which plotter draw.
			by default will be -1. \n
			0 : SciLab \n
			1 : GNU-Plot \n
			2 : MatLab \n
			3 : Maple  */

		virtual void check_prior(const Prob<eT>& pi) const {
			if(C.n_rows != pi.n_cols)
				throw std::runtime_error("invalid prior size");
		}
};


template<typename eT>
void LeakageMeasure<eT>::plot2d_vulnerability() {
	//correct size control
	if(C.n_rows != 2) {
		throw std::runtime_error("invalid size"); // X must be equal for both
	}

	//Should check for the specific classes if can implement it or not.
	//Shannon and Guessing have not vulnerability
	std::string className;
	className = this->class_name();
	if(className.compare("Shannon") == 0) {
		throw std::runtime_error("not implemented");
	}
	if(className.compare("Guessing") == 0) {
		throw std::runtime_error("not implemented");
	}
	switch(plotter_flag) {
	case -1:
		throw std::runtime_error("must choose engine");

	case 0: {
		std::ofstream outputFile;
		std::ifstream inputFile;
		std::string fileName;

		// creating the output file
		fileName = className + "_plot2d_vulnerability.sci";
		outputFile.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);

		// importing the vulnerability function to the file
		inputFile.open("../sci_files/2dvulnerability.sci", std::ifstream::in);
		std::string current_line;
		while(getline(inputFile, current_line)) {
			outputFile << current_line << std::endl;
		}
		inputFile.close();
		//writing channel C
		std::string channel_c = "C = [";
		//TODOchannel_c +=C.str;
		channel_c += "];";
		outputFile << channel_c << std::endl;

		//if it is GLeakage we should write Gain Function g
		if(className.compare("GLeakage") == 0) {
			std::string gain_g = "g = [";
			//TODOgain_g += g->str;
			gain_g += "];";
			outputFile << gain_g << std::endl;
		}

		// generating the plot
		std::string plotting_command = "[a,b]=size(pi);";
		outputFile << plotting_command << std::endl;

		plotting_command = "Res=zeros(b,a);";
		outputFile << plotting_command << std::endl;

		plotting_command = "for i=1:b";
		outputFile << plotting_command << std::endl;

		//choosing the correct vulnerability function
		if(className.compare("GLeakage") == 0) {
			plotting_command = "Res(i)=V(pi(i));";
			outputFile << plotting_command << std::endl;
		} else {
			//it is MinEntropy
			plotting_command = "Res(i)=Vm(pi(i));";
			outputFile << plotting_command << std::endl;
		}

		plotting_command = "end";
		outputFile << plotting_command << std::endl;

		plotting_command = "scf();";
		outputFile << plotting_command << std::endl;

		plotting_command = "plot(pi,Res);";
		outputFile << plotting_command << std::endl;
		// end of generating the plot

		outputFile.close();
	}
	break;

	case 1:
		throw std::runtime_error("not implemented");

	case 2:
		throw std::runtime_error("not implemented");

	case 3:
		throw std::runtime_error("not implemented");
	}
}

template<typename eT>
void LeakageMeasure<eT>::plot2d_cond_vulnerability() {
	//correct size control
	if(C.n_rows != 2) {
		throw std::runtime_error("invalid size"); // X must be equal for both
	}

	//Should check for the specific classes if can implement it or not.
	//Shannon and Guessing have not conditional vulnerability
	std::string className;
	className = this->class_name();
	if(className.compare("Shannon") == 0) {
		throw std::runtime_error("not implemented");
	}
	if(className.compare("Guessing") == 0) {
		throw std::runtime_error("not implemented");
	}
	switch(plotter_flag) {
	case -1:
		throw std::runtime_error("must choose engine");

	case 0: {
		std::ofstream outputFile;
		std::ifstream inputFile;
		std::string fileName;

		// creating the output file
		fileName = className + "_plot2d_cond_vulnerability.sci";
		outputFile.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);

		// importing the cond_vulnerability function to the file
		inputFile.open("../sci_files/2dcond_vulnerability.sci", std::ifstream::in);
		std::string current_line;
		while(getline(inputFile, current_line)) {
			outputFile << current_line << std::endl;
		}
		inputFile.close();
		//writing channel C
		std::string channel_c = "C = [";
		//TODOchannel_c +=C.str;
		channel_c += "];";
		outputFile << channel_c << std::endl;

		//if it is GLeakage we should write Gain Function g
		if(className.compare("GLeakage") == 0) {
			std::string gain_g = "g = [";
			//TODOgain_g += g->str;
			gain_g += "];";
			outputFile << gain_g << std::endl;
		}

		// generating the plot
		std::string plotting_command = "[a,b]=size(pi);";
		outputFile << plotting_command << std::endl;

		plotting_command = "Res=zeros(b,a);";
		outputFile << plotting_command << std::endl;

		plotting_command = "for i=1:b";
		outputFile << plotting_command << std::endl;

		//choosing the correct cond_vulnerability function
		if(className.compare("GLeakage") == 0) {
			plotting_command = "Res(i)=Vp(pi(i));";
			outputFile << plotting_command << std::endl;
		} else {
			//it is MinEntropy
			plotting_command = "Res(i)=Vpm(pi(i));";
			outputFile << plotting_command << std::endl;
		}

		plotting_command = "end";
		outputFile << plotting_command << std::endl;

		plotting_command = "scf();";
		outputFile << plotting_command << std::endl;

		plotting_command = "plot(pi,Res);";
		outputFile << plotting_command << std::endl;
		// end of generating the plot

		outputFile.close();
	}
	break;

	case 1:
	case 2:
	case 3:
		throw std::runtime_error("not implemented");
	}
}

template<typename eT>
void LeakageMeasure<eT>::plot2d_leakage() {
	//correct size control
	if(C.n_rows != 2) {
		throw std::runtime_error("invalid size"); // X must be equal for both
	}
	std::string className;
	className = this->class_name();
	switch(plotter_flag) {
	case -1:
		throw std::runtime_error("must choose engine");

	case 0: {
		std::ofstream outputFile;
		std::ifstream inputFile;
		std::string fileName;

		// creating the output file
		fileName = className + "_plot2d_leakage.sci";
		outputFile.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);

		// importing the leakage function to the file
		inputFile.open("../sci_files/2dleakage.sci", std::ifstream::in);
		std::string current_line;
		while(getline(inputFile, current_line)) {
			outputFile << current_line << std::endl;
		}
		inputFile.close();
		//writing channel C
		std::string channel_c = "C = [";
		//TODOchannel_c +=C.str;
		channel_c += "];";
		outputFile << channel_c << std::endl;

		//if it is GLeakage we should write Gain Function g
		if(className.compare("GLeakage") == 0) {
			std::string gain_g = "g = [";
			//TODOgain_g += g->str;
			gain_g += "];";
			outputFile << gain_g << std::endl;
		}

		// generating the plot
		std::string plotting_command = "[a,b]=size(pi);";
		outputFile << plotting_command << std::endl;

		plotting_command = "Res=zeros(b,a);";
		outputFile << plotting_command << std::endl;

		plotting_command = "for i=1:b";
		outputFile << plotting_command << std::endl;

		//choosing the correct leakage function
		if(className.compare("GLeakage") == 0) {
			plotting_command = "Res(i)=ML(pi(i));";
			outputFile << plotting_command << std::endl;
		} else {
			if(className.compare("MinEntropy") == 0) {
				//it is MinEntropy
				plotting_command = "Res(i)=MLm(pi(i));";
				outputFile << plotting_command << std::endl;
			} else {
				if(className.compare("Guessing") == 0) {
					//it is Guessing
					plotting_command = "Res(i)=MLg(pi(i));";
					outputFile << plotting_command << std::endl;
				} else {
					//it is Shannon
					plotting_command = "Res(i)=MLs(pi(i));";
					outputFile << plotting_command << std::endl;
				}
			}
		} // end if

		plotting_command = "end";
		outputFile << plotting_command << std::endl;

		plotting_command = "scf();";
		outputFile << plotting_command << std::endl;

		plotting_command = "plot(pi,Res);";
		outputFile << plotting_command << std::endl;
		// end of generating the plot

		outputFile.close();
	}
	break;

	case 1:
	case 2:
	case 3:
		throw std::runtime_error("not implemented");
	}
}

template<typename eT>
void LeakageMeasure<eT>::plot2d_entropy() {
	//correct size control
	if(C.n_rows != 2) {
		throw std::runtime_error("invalid size"); // X must be equal for both
	}
	std::string className;
	className = this->class_name();
	switch(plotter_flag) {
	case -1:
		throw std::runtime_error("must choose engine");

	case 0: {
		std::ofstream outputFile;
		std::ifstream inputFile;
		std::string fileName;

		// creating the output file
		fileName = className + "_plot2d_entropy.sci";
		outputFile.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);

		// importing the entropy function to the file
		inputFile.open("../sci_files/2dentropy.sci", std::ifstream::in);
		std::string current_line;
		while(getline(inputFile, current_line)) {
			outputFile << current_line << std::endl;
		}
		inputFile.close();
		//writing channel C
		std::string channel_c = "C = [";
		//TODOchannel_c +=C.str;
		channel_c += "];";
		outputFile << channel_c << std::endl;

		//if it is GLeakage we should write Gain Function g
		if(className.compare("GLeakage") == 0) {
			std::string gain_g = "g = [";
			//TODOgain_g += g->str;
			gain_g += "];";
			outputFile << gain_g << std::endl;
		}

		// generating the plot
		std::string plotting_command = "[a,b]=size(pi);";
		outputFile << plotting_command << std::endl;

		plotting_command = "Res=zeros(b,a);";
		outputFile << plotting_command << std::endl;

		plotting_command = "for i=1:b";
		outputFile << plotting_command << std::endl;

		//choosing the correct entropy function
		if(className.compare("GLeakage") == 0) {
			plotting_command = "Res(i)=E(pi(i));";
			outputFile << plotting_command << std::endl;
		} else {
			if(className.compare("MinEntropy") == 0) {
				//it is MinEntropy
				plotting_command = "Res(i)=Em(pi(i));";
				outputFile << plotting_command << std::endl;
			} else {
				if(className.compare("Guessing") == 0) {
					//it is Guessing
					plotting_command = "Res(i)=Eg(pi(i));";
					outputFile << plotting_command << std::endl;
				} else {
					//it is Shannon
					plotting_command = "Res(i)=Es(pi(i));";
					outputFile << plotting_command << std::endl;
				}
			}
		} // end if

		plotting_command = "end";
		outputFile << plotting_command << std::endl;

		plotting_command = "scf();";
		outputFile << plotting_command << std::endl;

		plotting_command = "plot(pi,Res);";
		outputFile << plotting_command << std::endl;
		// end of generating the plot

		outputFile.close();
	}
	break;

	case 1:
		//Its not implemented yet
		throw std::runtime_error("not implemented");
		break;

	case 2:
		//Its not implemented yet
		throw std::runtime_error("not implemented");
		break;

	case 3:
		//Its not implemented yet
		throw std::runtime_error("not implemented");
		break;
	}
}

template<typename eT>
void LeakageMeasure<eT>::plot2d_cond_entropy() {
	//correct size control
	if(C.n_rows != 2) {
		throw std::runtime_error("invalid size"); // X must be equal for both
	}
	std::string className;
	className = this->class_name();
	switch(plotter_flag) {
	case -1:
		throw std::runtime_error("must choose engine");

	case 0: {
		std::ofstream outputFile;
		std::ifstream inputFile;
		std::string fileName;

		// creating the output file
		fileName = className + "_plot2d_cond_entropy.sci";
		outputFile.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);

		// importing the cond_entropy function to the file
		inputFile.open("../sci_files/2dcond_entropy.sci", std::ifstream::in);
		std::string current_line;
		while(getline(inputFile, current_line)) {
			outputFile << current_line << std::endl;
		}
		inputFile.close();
		//writing channel C
		std::string channel_c = "C = [";
		//TODOchannel_c +=C.str;
		channel_c += "];";
		outputFile << channel_c << std::endl;

		//if it is GLeakage we should write Gain Function g
		if(className.compare("GLeakage") == 0) {
			std::string gain_g = "g = [";
			//TODOgain_g += g->str;
			gain_g += "];";
			outputFile << gain_g << std::endl;
		}

		// generating the plot
		std::string plotting_command = "[a,b]=size(pi);";
		outputFile << plotting_command << std::endl;

		plotting_command = "Res=zeros(b,a);";
		outputFile << plotting_command << std::endl;

		plotting_command = "for i=1:b";
		outputFile << plotting_command << std::endl;

		//choosing the correct cond_entropy function
		if(className.compare("GLeakage") == 0) {
			plotting_command = "Res(i)=Ep(pi(i));";
			outputFile << plotting_command << std::endl;
		} else {
			if(className.compare("MinEntropy") == 0) {
				//it is MinEntropy
				plotting_command = "Res(i)=Epm(pi(i));";
				outputFile << plotting_command << std::endl;
			} else {
				if(className.compare("Guessing") == 0) {
					//it is Guessing
					plotting_command = "Res(i)=Epg(pi(i));";
					outputFile << plotting_command << std::endl;
				} else {
					//it is Shannon
					plotting_command = "Res(i)=Eps(pi(i));";
					outputFile << plotting_command << std::endl;
				}
			}
		} // end if

		plotting_command = "end";
		outputFile << plotting_command << std::endl;

		plotting_command = "scf();";
		outputFile << plotting_command << std::endl;

		plotting_command = "plot(pi,Res);";
		outputFile << plotting_command << std::endl;
		// end of generating the plot

		outputFile.close();
	}
	break;

	case 1:
	case 2:
	case 3:
		throw std::runtime_error("not implemented");
	}
}

template<typename eT>
void LeakageMeasure<eT>::plot3d_vulnerability() {
	//correct size control
	if(C.n_rows != 3) {
		throw std::runtime_error("invalid size"); // X must be equal for both
	}

	std::string className;
	className = this->class_name();
	if(className.compare("Shannon") == 0) {
		throw std::runtime_error("not implemented");
	}
	if(className.compare("Guessing") == 0) {
		throw std::runtime_error("not implemented");
	}
	switch(plotter_flag) {
	case -1:
		throw std::runtime_error("must choose engine");

	case 0: {
		std::ofstream outputFile;
		std::ifstream inputFile;
		std::string fileName;

		// creating the output file
		fileName = className + "_plot3d_vulneability.sci";
		outputFile.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);

		// importing the vulnerability function to the file
		inputFile.open("../sci_files/3dvulneability.sci", std::ifstream::in);
		std::string current_line;
		while(getline(inputFile, current_line)) {
			outputFile << current_line << std::endl;
		}
		inputFile.close();
		//writing channel C
		std::string channel_c = "C = [";
		//TODOchannel_c +=C.str;
		channel_c += "];";
		outputFile << channel_c << std::endl;

		//if it is GLeakage we should write Gain Function g
		if(className.compare("GLeakage") == 0) {
			std::string gain_g = "g = [";
			//TODOgain_g += g->str;
			gain_g += "];";
			outputFile << gain_g << std::endl;
		}

		// generating the plot
		std::string plotting_command = "scf();";
		outputFile << plotting_command << std::endl;

		//choosing the correct vulnerability function
		if(className.compare("GLeakage") == 0) {
			plotting_command = "fplot3d(pi,pi,V,30,45,'X@Y@Z');";
			outputFile << plotting_command << std::endl;
		} else {
			//it is MinEntropy
			plotting_command = "fplot3d(pi,pi,Vm,30,45,'X@Y@Z');";
			outputFile << plotting_command << std::endl;
		} // end if
		// end of generating the plot

		outputFile.close();
	}
	break;

	case 1:
	case 2:
	case 3:
		throw std::runtime_error("not implemented");
	}
}

template<typename eT>
void LeakageMeasure<eT>::plot3d_cond_vulnerability() {
	//correct size control
	if(C.n_rows != 3) {
		throw std::runtime_error("invalid size"); // X must be equal for both
	}

	std::string className;
	className = this->class_name();
	if(className.compare("Shannon") == 0) {
		throw std::runtime_error("not implemented");
	}
	if(className.compare("Guessing") == 0) {
		throw std::runtime_error("not implemented");
	}
	switch(plotter_flag) {
	case -1:
		throw std::runtime_error("must choose engine");

	case 0: {
		std::ofstream outputFile;
		std::ifstream inputFile;
		std::string fileName;

		// creating the output file
		fileName = className + "_plot3d_cond_vulneability.sci";
		outputFile.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);

		// importing the cond_vulnerability function to the file
		inputFile.open("../sci_files/3dcond_vulneability.sci", std::ifstream::in);
		std::string current_line;
		while(getline(inputFile, current_line)) {
			outputFile << current_line << std::endl;
		}
		inputFile.close();
		//writing channel C
		std::string channel_c = "C = [";
		//TODOchannel_c +=C.str;
		channel_c += "];";
		outputFile << channel_c << std::endl;

		//if it is GLeakage we should write Gain Function g
		if(className.compare("GLeakage") == 0) {
			std::string gain_g = "g = [";
			//TODOgain_g += g->str;
			gain_g += "];";
			outputFile << gain_g << std::endl;
		}

		// generating the plot
		std::string plotting_command = "scf();";
		outputFile << plotting_command << std::endl;

		//choosing the correct cond_vulnerability function
		if(className.compare("GLeakage") == 0) {
			plotting_command = "fplot3d(pi,pi,Vp,30,45,'X@Y@Z');";
			outputFile << plotting_command << std::endl;
		} else {
			//it is MinEntropy
			plotting_command = "fplot3d(pi,pi,Vpm,30,45,'X@Y@Z');";
			outputFile << plotting_command << std::endl;
		} // end if
		// end of generating the plot

		outputFile.close();
	}
	break;

	case 1:
	case 2:
	case 3:
		throw std::runtime_error("not implemented");
	}
}

template<typename eT>
void LeakageMeasure<eT>::plot3d_leakage() {
	//correct size control
	if(C.n_rows != 3) {
		throw std::runtime_error("invalid size"); // X must be equal for both
	}

	std::string className;
	className = this->class_name();
	switch(plotter_flag) {
	case -1:
		throw std::runtime_error("must choose engine");

	case 0: {
		std::ofstream outputFile;
		std::ifstream inputFile;
		std::string fileName;

		// creating the output file
		fileName = className + "_plot3d_leakage.sci";
		outputFile.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);

		// importing the leakage function to the file
		inputFile.open("../sci_files/3dleakage.sci", std::ifstream::in);
		std::string current_line;
		while(getline(inputFile, current_line)) {
			outputFile << current_line << std::endl;
		}
		inputFile.close();
		//writing channel C
		std::string channel_c = "C = [";
		//TODOchannel_c +=C.str;
		channel_c += "];";
		outputFile << channel_c << std::endl;

		//if it is GLeakage we should write Gain Function g
		if(className.compare("GLeakage") == 0) {
			std::string gain_g = "g = [";
			//TODOgain_g += g->str;
			gain_g += "];";
			outputFile << gain_g << std::endl;
		}

		// generating the plot
		std::string plotting_command = "scf();";
		outputFile << plotting_command << std::endl;

		//choosing the correct leakage function
		if(className.compare("GLeakage") == 0) {
			plotting_command = "fplot3d(pi,pi,ML,30,45,'X@Y@Z');";
			outputFile << plotting_command << std::endl;
		} else {
			if(className.compare("MinEntropy") == 0) {
				//it is MinEntropy
				plotting_command = "fplot3d(pi,pi,MLm,30,45,'X@Y@Z');";
				outputFile << plotting_command << std::endl;
			} else {
				if(className.compare("Guessing") == 0) {
					//it is Guessing
					plotting_command = "fplot3d(pi,pi,MLg,30,45,'X@Y@Z');";
					outputFile << plotting_command << std::endl;
				} else {
					//it is Shannon
					plotting_command = "fplot3d(pi,pi,MLs,30,45,'X@Y@Z');";
					outputFile << plotting_command << std::endl;
				}
			}
		} // end if
		// end of generating the plot

		outputFile.close();
	}
	break;

	case 1:
	case 2:
	case 3:
		throw std::runtime_error("not implemented");
	}
}

template<typename eT>
void LeakageMeasure<eT>::plot3d_entropy() {
	//correct size control
	if(C.n_rows != 3) {
		throw std::runtime_error("invalid size"); // X must be equal for both
	}

	std::string className;
	className = this->class_name();
	switch(plotter_flag) {
	case -1:
		throw std::runtime_error("must choose engine");

	case 0: {
		std::ofstream outputFile;
		std::ifstream inputFile;
		std::string fileName;

		// creating the output file
		fileName = className + "_plot3d_entropy.sci";
		outputFile.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);

		// importing the entropy function to the file
		inputFile.open("../sci_files/3dentropy.sci", std::ifstream::in);
		std::string current_line;
		while(getline(inputFile, current_line)) {
			outputFile << current_line << std::endl;
		}
		inputFile.close();
		//writing channel C
		std::string channel_c = "C = [";
		//TODOchannel_c +=C.str;
		channel_c += "];";
		outputFile << channel_c << std::endl;

		//if it is GLeakage we should write Gain Function g
		if(className.compare("GLeakage") == 0) {
			std::string gain_g = "g = [";
			//TODOgain_g += g->str;
			gain_g += "];";
			outputFile << gain_g << std::endl;
		}

		// generating the plot
		std::string plotting_command = "scf();";
		outputFile << plotting_command << std::endl;

		//choosing the correct entropy function
		if(className.compare("GLeakage") == 0) {
			plotting_command = "fplot3d(pi,pi,E,30,45,'X@Y@Z');";
			outputFile << plotting_command << std::endl;
		} else {
			if(className.compare("MinEntropy") == 0) {
				//it is MinEntropy
				plotting_command = "fplot3d(pi,pi,Em,30,45,'X@Y@Z');";
				outputFile << plotting_command << std::endl;
			} else {
				if(className.compare("Guessing") == 0) {
					//it is Guessing
					plotting_command = "fplot3d(pi,pi,Eg,30,45,'X@Y@Z');";
					outputFile << plotting_command << std::endl;
				} else {
					//it is Shannon
					plotting_command = "fplot3d(pi,pi,Es,30,45,'X@Y@Z');";
					outputFile << plotting_command << std::endl;
				}
			}
		} // end if
		// end of generating the plot

		outputFile.close();
	}
	break;

	case 1:
	case 2:
	case 3:
		throw std::runtime_error("not implemented");
	}
}

template<typename eT>
void LeakageMeasure<eT>::plot3d_cond_entropy() {
	//correct size control
	if(C.n_rows != 3) {
		throw std::runtime_error("invalid size"); // X must be equal for both
	}

	std::string className;
	className = this->class_name();
	switch(plotter_flag) {
	case -1:
		throw std::runtime_error("must choose engine");

	case 0: {
		std::ofstream outputFile;
		std::ifstream inputFile;
		std::string fileName;

		// creating the output file
		fileName = className + "_plot3d_cond_entropy.sci";
		outputFile.open(fileName.c_str(), std::ofstream::out | std::ofstream::app);

		// importing the cond_entropy function to the file
		inputFile.open("../sci_files/3dcond_entropy.sci", std::ifstream::in);
		std::string current_line;
		while(getline(inputFile, current_line)) {
			outputFile << current_line << std::endl;
		}
		inputFile.close();
		//writing channel C
		std::string channel_c = "C = [";
		//TODOchannel_c +=C.str;
		channel_c += "];";
		outputFile << channel_c << std::endl;

		//if it is GLeakage we should write Gain Function g
		if(className.compare("GLeakage") == 0) {
			std::string gain_g = "g = [";
			//TODOgain_g += g->str;
			gain_g += "];";
			outputFile << gain_g << std::endl;
		}

		// generating the plot
		std::string plotting_command = "scf();";
		outputFile << plotting_command << std::endl;

		//choosing the correct cond_entropy function
		if(className.compare("GLeakage") == 0) {
			plotting_command = "fplot3d(pi,pi,Ep,30,45,'X@Y@Z');";
			outputFile << plotting_command << std::endl;
		} else {
			if(className.compare("MinEntropy") == 0) {
				//it is MinEntropy
				plotting_command = "fplot3d(pi,pi,Epm,30,45,'X@Y@Z');";
				outputFile << plotting_command << std::endl;
			} else {
				if(className.compare("Guessing") == 0) {
					//it is Guessing
					plotting_command = "fplot3d(pi,pi,Epg,30,45,'X@Y@Z');";
					outputFile << plotting_command << std::endl;
				} else {
					//it is Shannon
					plotting_command = "fplot3d(pi,pi,Eps,30,45,'X@Y@Z');";
					outputFile << plotting_command << std::endl;
				}
			}
		} // end if
		// end of generating the plot

		outputFile.close();
	}
	break;

	case 1:
	case 2:
	case 3:
		throw std::runtime_error("not implemented");
	}
}

template<typename eT>
void LeakageMeasure<eT>::change_to_scilab() {
	switch(plotter_flag) {
	case -1:
		plotter_flag = 0;
		//open scilab
		break;

	case 0: //do nothing
		break;

	case 1: //close gnuplot
		plotter_flag = -1;
		change_to_scilab();
		break;

	case 2: //close matlab
		plotter_flag = -1;
		change_to_scilab();
		break;

	case 3: //close maple
		plotter_flag = -1;
		change_to_scilab();
		break;
	}
}

template<typename eT>
void LeakageMeasure<eT>::change_to_gnuplot() {
	throw std::runtime_error("not implemented");
}

template<typename eT>
void LeakageMeasure<eT>::change_to_matlab() {
	throw std::runtime_error("not implemented");
}

template<typename eT>
void LeakageMeasure<eT>::change_to_maple() {
	throw std::runtime_error("not implemented");
}
