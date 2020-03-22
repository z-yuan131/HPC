#include "Poisson.h"
using namespace std;
#include <iostream>
#include <string>
#include <fstream> //write to files
#include <algorithm>    // std::max
#include <math.h>       /* fabs sqrt*/
#include <boost/program_options.hpp>
#include <mpi.h>
// namespace po = boost::program_options;

// #include "cblas.h"

#include <cstdlib>

#define F77Name(x) x##_
extern "C"{
    //LAPACK routine for solving systems of linear equations
    void F77Name(dgbsv)(const int& n, const int&kl, const int& ku,
                        const int& nrhs, const double * AB,        //B stands for banded
                        const int& ldab, int * ipiv, double * B,    //ipiv: vector for pivot
                        const int& ldb, int& info);
}
void Poisson::set(LidDrivenCavity* Lid){
  v_in  = Lid->v_in;
  s_in  = Lid->s_in;
  s_bcB = Lid->s_bcB;
  s_bcT = Lid->s_bcT;
  s_bcL = Lid->s_bcL;
  s_bcR = Lid->s_bcR;
  A     = Lid->A;
  n     = Lid->dNx*Lid->dNy;    //matrix A size
  kl    = Lid->dNy;    //Lower diagonal bandwidth
  ku    = Lid->dNy;    //Upper diagonal bandwidth
  nrhs  = 1;     //number of right hand side vectors
  ldab  = 1 + 2*kl + ku;   //Number of rows in compressed matrix for Lapack
  info = 1;
}


void Poisson::PoissonSolver(LidDrivenCavity* Lid){


     double* AB = new double[ldab*n]; //for Lapack will overwrite matrix A;
     int* ipiv = new int[n];
     double* temp = new double[n];      //for lapack will overwrite v with s;


     for (int i = 0 ; i < ldab*n; i++){
       AB[i] = A[i];
     }

     for (int i = 0 ; i < n ; i++){
       temp[i] = v_in[i];
     }

     //apply boundary conditions
     for (int j = 0 ; j < Lid->dNy ; j++){
       for (int i = 0 ; i < Lid->dNx ; i++){
         if (j == 0){         //bottom bc
           temp[j + Lid->dNy*i] += s_bcB[i]/Lid->dy/Lid->dy;
           if (i == 0){
             temp[j + Lid->dNy*i] += s_bcL[j]/Lid->dx/Lid->dx;
           }
           else if(i == Lid->dNx - 1){
             temp[j + Lid->dNy*i] += s_bcR[j]/Lid->dx/Lid->dx;
           }


         }
         else if(j == Lid->dNy - 1){      //top bc
           temp[j + Lid->dNy*i] += s_bcT[i]/Lid->dy/Lid->dy;

           if (i == 0){
             temp[j + Lid->dNy*i] += s_bcL[j]/Lid->dx/Lid->dx;
           }
           else if(i == Lid->dNx - 1){
             temp[j + Lid->dNy*i] += s_bcR[j]/Lid->dx/Lid->dx;
           }

         }
         else if(i == 0){    //left bc
           temp[j + Lid->dNy*i] += s_bcL[j]/Lid->dx/Lid->dx;
         }
         else if(i == Lid->dNx - 1){      //right bc
           temp[j + Lid->dNy*i] += s_bcR[j]/Lid->dx/Lid->dx;
         }
       }
     }


       // for (int i = 0 ; i < dNy ; i++){
       //   cout <<" R_s "<< s_bcR[i] << "  ";
       // }cout << endl;






     // double* error_temp = new double[dNx*dNy];
     // for (int i = 0 ; i < n ; i++){
     //   error_temp[i] = s_in[i];
     // }
     // s_in_error = error_temp;

     // cout << "Poisson solver:AB" << endl;
     // for (int j = 0 ; j < ldab ; j++){
     //   for (int i = 0; i < n ; i++){
     //     cout << AB[j + i*ldab] << "    ";
     //   }
     //   cout << "\n";
     // }
     // cout << endl;
     //
     // cout << "Poisson solver:vorticity apply bc" << endl;
     // for (int j = 0 ; j < dNy ; j++){
     //   for (int i = 0; i < dNx ; i++){
     //     cout << v_in[j + i*dNy] << "    ";
     //   }
     //   cout << "\n";
     // }
     // cout << "interior vorticity before" << endl;
     //   for (int i = 0; i < n ; i++){
     //     cout << temp[i] << "    ";
     //   }cout << endl;


     // s_in_error = s_in;

     F77Name(dgbsv)(n, kl, ku, nrhs, AB, ldab, ipiv, temp, n, info); //temp is input and output
     // cout << "info = " << info << endl;
     // cout << "interior vorticity after" << endl;
     // for (int j = 0 ; j < dNy ; j++){
     //   for (int i = 0; i < dNx ; i++){
     //     cout << temp[j + i*dNy] << "    ";
     //   }
     //   cout << "\n";
     // }
     for (int i = 0 ; i < n ; i++){
       s_in[i] = temp[i];
     }
     // s_in = temp;
     // for (int j = 0 ; j < dNy ; j++){
     //   for (int i = 0; i < dNx ; i++){
     //     cout << s_in[j + i*dNy] << "    ";
     //   }
     //   cout << "\n";
     // }

     delete[] temp;   //bug
     delete[] AB;
     delete[] ipiv;


     // cout << "Poisson solver:streamFunction" << endl;
     // // cout << "streamFunction" << s_in[10] << endl;
     // for (int j = 0 ; j < dNy ; j++){
     //   for (int i = 0; i < dNx ; i++){
     //     cout << s_in[j + i*dNy] << "    ";
     //   }
     //   cout << "\n";
     // }
     // cout << endl;

 }
