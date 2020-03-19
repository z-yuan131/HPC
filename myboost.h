#pragma once
#include <string>

class myboost{
public:
  //read input function
  void myboostfun(int argc, char **argv);

  //return functions
  int returnNx();
  int returnNy();
  int returnPx();
  int returnPy();

  double returnLx();
  double returnLy();
  double returndt();
  double returnT();
  double returnRe();


private:
  int Nx;
  int Ny;
  int Px;
  int Py;
  double Lx;
  double Ly;
  double dt;
  double T ;
  double Re;
};
