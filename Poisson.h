#pragma once
#include "LidDrivenCavity.h"
#include <string>

class LidDrivenCavity;

class Poisson
{


public:
  // LidDrivenCavity* Lid;
  void set(LidDrivenCavity* Lid);
  void PoissonSolver(LidDrivenCavity* Lid);

private:
  double* v_in = nullptr;      //array of interior vorticity
  double* s_in = nullptr;      //array of interior stream function
  double* s_bcT = nullptr;      //array of boundary stream function
  double* s_bcB = nullptr;
  double* s_bcL = nullptr;      //array of boundary stream function
  double* s_bcR = nullptr;
  double* A = nullptr;      //matrix(array) of linear system Ax = y

};
