#include <iostream>
#include <armadillo>
using namespace arma;

int main()
{
	
    std::cout << "Using Armadillo 1" << std::endl;
    char* channel_elements =(char *)malloc((strlen("1 0 0")+1)*sizeof(char));
    std::cout << "Using Armadillo 2" << std::endl;
    strcpy(channel_elements, "1 0 0");
    std::cout << "Using Armadillo 3" << std::endl;
    //vec v=vec(channel_elements);
    mat v;
    std::cout << "Using Armadillo 4" << std::endl;
    return 0;
}
