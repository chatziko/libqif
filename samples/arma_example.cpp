#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main(int argc, char** argv){

  char* channel_elements =(char *)malloc((strlen("1 0 0")+1)*sizeof(char));
  strcpy(channel_elements, "1 0 0");	

  mat A = randu<mat>(4,5);
  mat B = randu<mat>(4,5);
  
  cout << A*B.t() << endl;
  
  return 0;
}
