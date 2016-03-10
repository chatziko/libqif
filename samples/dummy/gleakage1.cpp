#include <iostream>
#include <string>
#include "GLeakage.h"

void gleakage1() {
	std::cout << "Using LIBQIF Library looking for properties" << std::endl;

	std::string random = "0.3 0.7; 0.7 0.3; 0.3 0.7";
	std::string balanced = "0.25 0.5 0.25; 0 1 0; 0.5 0 0.5";
	std::string metrics = "1 0.5 0; 0.5 1 0.5; 0 0.5 1";
	std::string id = "1 0 0; 0 1 0; 0 0 1";
	std::string k_tries = "1 1 0; 1 0 1; 0 1 1";

	chan C_rand = chan(random);
	chan C_balanced = chan(balanced);
	chan C_id = chan(id);

	mat g_id = mat(id);
	mat g_metrics = mat(metrics);
	mat g_k_tries = mat(k_tries);

	//GLeakage<double>
	GLeakage<double> gl1 = GLeakage<double>(C_rand, g_id);
	GLeakage<double> gl2 = GLeakage<double>(C_balanced, g_id);
	GLeakage<double> gl3 = GLeakage<double>(C_id, g_id);

	GLeakage<double> gl4 = GLeakage<double>(C_rand, g_metrics);
	GLeakage<double> gl5 = GLeakage<double>(C_balanced, g_metrics);
	GLeakage<double> gl6 = GLeakage<double>(C_id, g_metrics);

	GLeakage<double> gl7 = GLeakage<double>(C_rand, g_k_tries);
	GLeakage<double> gl8 = GLeakage<double>(C_balanced, g_k_tries);
	GLeakage<double> gl9 = GLeakage<double>(C_id, g_k_tries);

	gl1.change_to_scilab();
	gl2.change_to_scilab();
	gl3.change_to_scilab();

	gl4.change_to_scilab();

	gl5.change_to_scilab();
	gl6.change_to_scilab();
	gl7.change_to_scilab();
	gl8.change_to_scilab();
	gl9.change_to_scilab();

	//ploting

	gl1.plot3d_leakage();
	gl2.plot3d_leakage();
	gl3.plot3d_leakage();

	gl4.plot3d_leakage(); //<-----

	gl5.plot3d_leakage();
	gl6.plot3d_leakage();
	gl7.plot3d_leakage();
	gl8.plot3d_leakage();
	gl9.plot3d_leakage();

}
