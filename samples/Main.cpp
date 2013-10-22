#include <iostream>
#include <string>

#include "Shannon.h"
#include "MinEntropy.h"
#include "Guessing.h"
int main()
{
	
    std::cout << "Using QIF Library Example" << std::endl;
   
    
    // NOT SOPPORTED FOR THE MOMENT
    /*double channel_elements[3][3]={{1,0,0},{0,1,0},{0,0,1}};
    (double[3])* point= {{1,0,0},{0,1,0},{0,0,1}};
    */
    //IN JUST ONE LINE: Channel C= Channel({{1,0,0},{0,1,0},{0,0,1}});
    //IN JUST ONE LINE: Gain_function g=Gain_function("1 0 0; 0 1 0; 0 0 1");
    
    //Creating the channel matrix
    char* channel_elements =(char *)malloc((strlen("1 0 0; 0 1 0; 0 0 1")+1)*sizeof(char));
    strcpy(channel_elements, "1 0 0; 0 1 0; 0 0 1");

    
    Channel C= Channel((char*)"1 0 0; 0 1 0; 0 0 1");

    std::string pepe = "1 0 0; 0 1 0; 0 0 1"; 
    Channel C2= Channel(pepe);

    

    //Creating the gain function matrix
    //char* gain_function_elements =(char *)malloc((strlen("1 0 0; 0 1 0; 0 0 1")+1)*sizeof(char));
    //strcpy(gain_function_elements, "1 0 0; 0 1 0; 0 0 1");
    //Gain g=Gain(gain_function_elements);
    
    //Creating the probability distribution
    //double vector_elements[3]={1/3,1/3,1/3};
    char* vector_elements =(char *)malloc((strlen("1/3 1/3 1/3")+1)*sizeof(char));
    strcpy(vector_elements, "1/3 1/3 1/3");
    Prob p1= Prob(vector_elements);
    std::cout << "Calculating the GLeakage" << std::endl;
    //Calculating the GLeakage
    //GLeakage gl= GLeakage(C,g);
    Shannon sl = Shannon(C);
    MinEntropy ml = MinEntropy(C);
    Guessing gul = Guessing(C);
    std::cout << "Calculating ends" << std::endl;

    //double Lg=gl.leakage(p1);
    //double Ls=sl.leakage(p1);
    
    //std::cout << "Lg " << Lg << std::endl;
    //std::cout << "Ls " << Ls << std::endl;
     
    //Using the Plotter 
    //gl.plot3d_vulnerability();
    sl.change_to_scilab();
    sl.plot3d_leakage();
    
    
    //Using Linear Programming
    //LinearProgram lp;
    //vec v= lp.solve("1 1 1; 10 4 5; 2 2 6","","10 6 4","-100; -600; -300");
}
