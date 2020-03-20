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
  int     n;    //matrix A size
  int     kl;    //Lower diagonal bandwidth
  int     ku;    //Upper diagonal bandwidth
  int     nrhs;     //number of right hand side vectors
  int     ldab;   //Number of rows in compressed matrix for Lapack
  int     info;

};
