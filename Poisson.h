#pragma once
#include "LidDrivenCavity.h"      //this class is friend of LidDrivenCavity
#include <string>

class LidDrivenCavity;

class Poisson
{


public:
  void set(LidDrivenCavity* Lid);           //get information from class LidDrivenCavity
  void LUfact();                            //pre LU
  void PoissonSolver(LidDrivenCavity* Lid);

private:
  double* v_in = nullptr;      //array of interior vorticity
  double* s_in = nullptr;      //array of interior stream function
  double* s_bcT = nullptr;      //array of boundary stream function
  double* s_bcB = nullptr;
  double* s_bcL = nullptr;      //array of boundary stream function
  double* s_bcR = nullptr;
  double* A = nullptr;      //matrix(array) of linear system Ax = y
  double* AB = nullptr;
  double* temp = nullptr;
  int   * ipiv = nullptr;

  int     n;    //matrix A size
  int     kl;    //Lower diagonal bandwidth
  int     ku;    //Upper diagonal bandwidth
  int     nrhs;     //number of right hand side vectors
  int     ldab;   //Number of rows in compressed matrix for Lapack
  int     info;

};
