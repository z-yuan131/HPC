 #include "LidDrivenCavity.h"
 // #include "PoissonSolver.h"
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
   void F77Name(dgesv)(const int& n, const int& nrhs, const double * A,
                       const int& lda, int * ipiv, double * B,    //ipiv: vector for pivot
                       const int& ldb, int& info);
   //LAPACK routine for solving systems of linear equations
   void F77Name(dgbsv)(const int& n, const int&kl, const int& ku,
                       const int& nrhs, const double * AB,        //B stands for banded
                       const int& ldab, int * ipiv, double * B,    //ipiv: vector for pivot
                       const int& ldb, int& info);
   //BLAS routine for performing matrix-vector multiplication
   void F77Name(dgemv)(const char& trans, const int& m,
                       const int& n,      const double& alpha,
                       const double* A,   const int& lda,
                       const double* x,   const int& incx,
                       const double& beta,double* y,
                       const int& incy);
   //BLAS routine for performing banded matrix-vector multiplication
   void F77Name(dgbmv)(const char& trans, const int& m,
                       const int& n,      const int& kl,
                       const int& ku,     const double& alpha,
                       const double* A,   const int& lda,
                       const double* x,   const int& incx,
                       const double& beta,double* y,
                       const int& incy);
   //BLAS routine for performing symmtric banded matrix-vector multiplication
   void F77Name(dsbmv)(const char& UPLO,  //UPLO is for upper(U) or lower(L) part of A
                       const int& n,      const int& k,
                       const double& alpha,
                       const double* A,   const int& lda,
                       const double* x,   const int& incx,
                       const double& beta,double* y,
                       const int& incy);
 }


 LidDrivenCavity::LidDrivenCavity(double Lx, double Ly, int Nx, int Ny, int Px, int Py, double dt, double T, double Re)
 {
   // cout << Lx << Ly << Nx << Ny << Px << Py << dt << T <<Re << endl;
 }

 LidDrivenCavity::LidDrivenCavity()
 {
   cout << "it is default constructor" << endl;
 }

 LidDrivenCavity::~LidDrivenCavity()
 {
   cout << "eliminate default constructor" << endl << endl;
 }

 void LidDrivenCavity::SetDomainSize(double xlen, double ylen)
 {
     Lx = xlen;
     Ly = ylen;
     // cout << xlen << endl;
 }

 void LidDrivenCavity::SetGridSize(int nx, int ny)
 {
   Nx = nx;
   Ny = ny;
 }

 void LidDrivenCavity::SetTimeStep(double deltat)
 {
   dt = deltat;
 }

 void LidDrivenCavity::SetFinalTime(double finalt)
 {
   T = finalt;
 }

 void LidDrivenCavity::SetReynoldsNumber(double re)
 {
   Re = re;
 }

 void LidDrivenCavity::Setmeshsize(double lx, double ly, int nx, int ny)
 {
   dx = lx/(nx-1);
   dy = ly/(ny-1);
 }

 void LidDrivenCavity::SetPxPy(int px, int py){
   Px = px;
   Py = py;
 }

 void LidDrivenCavity::SetSubdomainGrids(int Nx, int Ny){
     dNx = (Nx-2) / Px;
     dNy = (Ny-2) / Py;
   // cout << "dNx = " << dNx <<endl;
   // cout << "dNy = " << dNy <<endl;
 }

 void LidDrivenCavity::Initialise()
 {

     double* phi = new double[dNx*dNy];
     double* omega = new double[dNx*dNy];
     for (int j = 0 ; j < dNy ; j++){
       for (int i = 0 ; i < dNx ; i++){
         phi[j + dNy*i] = 0;
         omega[j + dNy*i] = 0;
       }
     }
     double* phi_bcL = new double[dNy];
     double* phi_bcR = new double[dNy];
     double* phi_bcB = new double[dNx];
     double* phi_bcT = new double[dNx];

     for (int i = 0 ; i < dNx ; i++){
       phi_bcT[i] = 0;           //1-2lines: top-bottom;
       phi_bcB[i] = 0;
     }
     for (int j = 0 ; j < dNy ; j++){
       phi_bcL[j] = 0;           //1-2lines: left-right;
       phi_bcR[j] = 0;
     }


     // cout << phi[dNy-3 + dNy*(dNx-2)] << endl;
     s_in = phi;
     v_in = omega;
     s_bcL = phi_bcL;
     s_bcR = phi_bcR;
     s_bcB = phi_bcB;
     s_bcT = phi_bcT;

     // cout << "interior streamFunction initial" << endl;
     // for (int j = 0 ; j < dNy ; j++){
     //   for (int i = 0; i < dNx ; i++){
     //     cout << s_in[j + i*dNy] << "    ";
     //   }
     //   cout << "\n";
     // }


       // for (int i = 0 ; i < dNx ; i++){
       //   cout << s_bcT[i] << "  ";
       // }cout << endl;
       // for (int i = 0 ; i < dNx ; i++){
       //     cout << s_bcB[i] << "  ";
       //   }cout << endl;
       // for (int i = 0 ; i < dNy ; i++){
       //     cout << s_bcR[i] << "  ";
       //   }cout << endl;
       // for (int i = 0 ; i < dNy ; i++){
       //     cout << s_bcL[i] << "  ";
       //   }cout << endl;

     // delete[] phi;
     // delete[] omega;
     // delete[] ohi_bc;


 }

 /*
 void LidDrivenCavity::Initialise_top_boundary(int sizeNx, int sizeNy, double deltay)
 {
     double* phi = new double[sizeNx*sizeNy];
     double* omega = new double[sizeNx*sizeNy];
     for (int i = 0 ; i < sizeNx ; i++){
       for (int j = 0 ; j < sizeNy ; j++){
         phi[j + sizeNx*i] = 0;
         if (i == 0){
           omega[j + sizeNx*i] = 1.0/deltay;}
         else {
           omega[j + sizeNx*i] = 0;}
       }
     }
     // cout << phi[2] << endl;
     s = phi;
     v = omega;
     delete[] phi;
     delete[] omega;
 }*/

 void LidDrivenCavity::Integrate()
 {
 }


 void LidDrivenCavity::CalculateVorticityBC(int TopBC)
 {//from equations (6)(7)(8)(9)
    // int dNy = Ny - 2;
    // int dNx = Nx - 2;
    double* tempL = new double[dNy];
    double* tempR = new double[dNy];
    double* tempB = new double[dNx];
    double* tempT = new double[dNx];
    for (int j = 0 ; j < dNy ; j++){
      for (int i = 0 ; i < dNx ; i++){
        if (i == 0){                  //left
            tempL[j] = (s_bcL[j] - s_in[j + dNy*0])*2.0/dx/dx;
        }
        if (i == dNx-1){               //right
            tempR[j] = (s_bcR[j] - s_in[j + dNy*(dNx-1)])*2.0/dx/dx;
        }
        if (j == 0){                  //bottom
            tempB[i] = (s_bcB[i] - s_in[0 + dNy*i])*2.0/dy/dy;
        }
        if (j == dNy-1){               //top
            tempT[i] = (s_bcT[i] - s_in[dNy-1 + dNy*i])*2.0/dy/dy - TopBC*2.0*1.0/dy; //U = 1.0
        }
      }
    }
  v_bcL = tempL;
  v_bcR = tempR;
  v_bcB = tempB;
  v_bcT = tempT;

        /*for (int i = 0 ; i < dNx ; i++){
            cout << v_bcT[i] << "T  ";
          }cout << endl;
        for (int i = 0 ; i < dNx ; i++){
            cout << v_bcB[i] << "B  ";
          }cout << endl;
        for (int i = 0 ; i < dNy ; i++){
            cout << v_bcR[i] << "R  ";
          }cout << endl;
        for (int i = 0 ; i < dNy ; i++){
            cout << v_bcL[i] << "L  ";
          }cout << endl;*/

  // delete[] temp;
 }


 void LidDrivenCavity::BuildMatrixA_B_C()    //dNx: points of interior domain.
 {
     // int dNx = Nx-2;
     // int dNy = Ny-2;
     int n = dNx*dNy;    //matrix A size
     int kl = dNy;    //Lower diagonal bandwidth
     int ku = dNy;    //Upper diagonal bandwidth
     int nrhs = 1;     //number of right hand side vectors
     int ldab = 1 + 2*kl + ku;   //Number of rows in compressed matrix for Lapack
     int bdab = 1 + kl + ku;   //Number of rows in compressed matrix for Blas
     int ldb = n;      //size of RHS vector
     double* A_vv   = new double[bdab*n];
     double* A_temp = new double[ldab*n]; //for Lapack will rewrite A matrix
     double* B_temp = new double[3*n]; //B matrix uses blas routine
     double* C_temp = new double[bdab*n]; //C matrix uses blas routine

     //populate matrix A, C here:
     for (int i = 0 ; i < n ; i++){
       for (int j = 0; j < ldab ; j++){
         // for lapack calculation
         if (j == kl+ku){        //for diagonal
           A_temp[j + i*ldab] = 2.0/(dx*dx) + 2.0/(dy*dy);
         }
         else if (j == kl+ku-1 && (i%dNy)){   //for first upper diag
           A_temp[j + i*ldab] = -1.0/(dy*dy);
         }
         else if (j == kl && i >= dNy){   //for second upper diag
           A_temp[j + i*ldab] = -1.0/(dx*dx);
         }
         else if (j == kl+ku+1 && ((i+1)%dNy)){   //for first lower diag
           A_temp[j + i*ldab] = -1.0/(dy*dy);
         }
         else if (j == kl+ku+kl && i < n - dNy){   //for second lower diag
           A_temp[j + i*ldab] = -1.0/(dx*dx);
         }
         else {
           A_temp[j + i*ldab] = 0.0;
         }
         // for blas calculation
         if (j < bdab){
           if (j == ku){        //for diagonal
             A_vv[j + i*bdab] = 2.0/dx/dx + 2.0/dy/dy;
           }
           else if (j == ku-1 && (i%dNy)){   //for first upper diag
             A_vv[j + i*bdab] = -1.0/dy/dy;
           }
           else if (j == 0 && i >= dNy){   //for second upper diag
             A_vv[j + i*bdab] = -1.0/dx/dx;
             C_temp[j + i*bdab] = 1.0/2.0/dx;
           }
           else if (j == ku+1 && ((i+1)%dNy)){   //for first lower diag
             A_vv[j + i*bdab] = -1.0/dy/dy;
           }
           else if (j == ku+kl && i < n - dNy){   //for second lower diag
             A_vv[j + i*bdab] = -1.0/dx/dx;
             C_temp[j + i*bdab] = -1.0/2.0/dx;
           }
           else {
             A_vv[j + i*bdab] = 0.0;
             C_temp[j + i*bdab] = 0.0;
           }
         }
       }
     }
     A = A_temp;
     A_v = A_vv;
     // delete[] A_vv;
     // delete[] A_temp;

     //populate matrix B here:
     for (int i = 0 ; i < n ; i++){
       for (int j = 0; j < 3 ; j++){
         if (j == 0 && (i%dNy)){
           B_temp[j + i*3] = 1.0/2.0/dy;
         }
         else if (j == 2 && ((i+1)%dNy)){
           B_temp[j + i*3] = -1.0/2.0/dy;
         }
         else {
           B_temp[j + i*3] = 0.0;
         }
       }
     }

     // for (int i = 0 ; i < n ; i++){
     //   for (int j = 0; j < 1+ku ; j++){
     //     if (j == ku){        //for diagonal
     //       A_vv[j + i*bdab] = 2.0/dx/dx + 2.0/dy/dy;
     //     }
     //     else if (j == ku-1 && (i%dNy)){   //for first upper diag
     //       A_vv[j + i*bdab] = -1.0/dx/dx;
     //     }
     //     else if (j == 0 && i >= dNy){   //for second upper diag
     //       A_vv[j + i*bdab] = -1.0/dy/dy;
     //       C_temp[j + i*bdab] = 1.0/2.0/dx;
     //     }
     //     else {
     //       A_vv[j + i*bdab] = 0.0;
     //       C_temp[j + i*bdab] = 0.0;
     //     }
     //   }
     // }


     C = C_temp;
     B = B_temp;

     /*
     cout << "for LAPACK" << endl;
     for (int j = 0 ; j < ldab ; j++){
       for (int i = 0; i < n ; i++){
         cout << A[j + i*ldab] << "    ";
       }
       cout << "\n";
     }
     cout << "for BLAS" << endl;
     for (int j = 0 ; j < bdab ; j++){
       for (int i = 0; i < n ; i++){
         cout << A_v[j + i*bdab] << "    ";
       }
       cout << "\n";
     }


     cout << "for B mAtrix" << endl;
     for (int j = 0 ; j < 3 ; j++){
       for (int i = 0; i < n ; i++){
         cout << B[j + i*3] << "    ";
       }
       cout << "\n";
     }

     cout << "for C Matrix" << endl;
     for (int j = 0 ; j < bdab ; j++){
       for (int i = 0; i < n ; i++){
         cout << C[j + i*bdab] << "    ";
       }
       cout << "\n";
     }*/
 }




 void LidDrivenCavity::CalculateInteriorVorticityAtTimet()
 {//from equation (10)
     int kl = dNy;    //Lower diagonal bandwidth
     int ku = dNy;    //Upper diagonal bandwidth
     int n = dNx*dNy;    //matrix A size
     int bdab = 1 + kl + ku;   //Number of rows in compressed matrix for Blas

     //v = A_v*s;  vorticity = coeff-matrix*streamFunctiom ,matrix is banded
     // cout << "for BLAS" << endl;
     // for (int j = 0 ; j < bdab ; j++){
     //   for (int i = 0; i < n ; i++){
     //     cout << A_v[j + i*bdab] << "    ";
     //   }
     //   cout << "\n";
     // }
     // for (int i = 0 ; i < n ; i++){
     //   s_in[i] = 1.0;
     // }
     // s_in[6] = 0.1;
     // cout << "interior vorticity before" << endl;
     //   for (int i = 0; i < n ; i++){
     //     cout << s_in[i] << "    ";
     //   }cout << endl;

     F77Name(dgbmv)('N', n, n, kl, ku, 1.0, A_v, bdab, s_in, 1, 0.0, v_in, 1);
     //
     // cout << "interior vorticity after" << endl;
     // for (int j = 0 ; j < dNy ; j++){
     //   for (int i = 0; i < dNx ; i++){
     //     cout << v_in[j + i*dNy] << "    ";
     //   }
     //   cout << "\n";
     // }
     // F77Name(dsbmv)('U', n, ku, 1.0, A_v, 1+ku, s_in, 1, 0.0, v_in, 1);
     // cout << "interior vorticity after" << endl;
     // for (int j = 0 ; j < dNy ; j++){
     //   for (int i = 0; i < dNx ; i++){
     //     cout << v_in[j + i*dNy] << "    ";
     //   }
     //   cout << "\n";
     // }


     //apply boundary conditions
     for (int j = 0 ; j < dNy ; j++){
       for (int i = 0 ; i < dNx ; i++){
         if (j == 0){         //bottom bc
           if(i == 0){             //conner
             v_in[j + dNy*i] -= (s_bcL[j]/dx/dx + s_bcB[i]/dy/dy);
           }
           else if(i == dNx - 1){  //conner
             v_in[j + dNy*i] -= (s_bcR[j]/dx/dx + s_bcB[i]/dy/dy);
           }
           else{
             v_in[j + dNy*i] -= s_bcB[i]/dy/dy;
           }
         }
         else if(j == dNy - 1){      //top bc
           if(i == 0){      //conner
             v_in[j + dNy*i] -= (s_bcL[j]/dx/dx + s_bcT[i]/dy/dy);
           }
           else if(i == dNx - 1){  //conner
             v_in[j + dNy*i] -= (s_bcR[j]/dx/dx + s_bcT[i]/dy/dy);
           }
           else{
             v_in[j + dNy*i] -= s_bcT[i]/dy/dy;
           }
         }
         else if(i == 0){    //left bc
           v_in[j + dNy*i] -= s_bcL[j]/dx/dx;
         }
         else if(i == dNx - 1){      //right bc
           v_in[j + dNy*i] -= s_bcR[j]/dx/dx;
         }
       }
     }


     /*
     cout << "interior streamFunction at time t" << endl;
     for (int j = 0 ; j < dNy ; j++){
       for (int i = 0; i < dNx ; i++){
         cout << s_in[j + i*dNy] << "    ";
       }
       cout << "\n";
     }

     cout << "interior vorticity at time t" << endl;
     for (int j = 0 ; j < dNy ; j++){
       for (int i = 0; i < dNx ; i++){
         cout << v_in[j + i*dNy] << "    ";
       }
       cout << "\n";
     }*/


 }


 void LidDrivenCavity::TimeAdvance(){
     int n = dNx*dNy;    //matrix A size
     int kl = dNy;    //Lower diagonal bandwidth
     int ku = dNy;    //Upper diagonal bandwidth
     int klB = 1;    //Lower diagonal bandwidth
     int kuB = 1;    //Upper diagonal bandwidth
     int nrhs = 1;     //number of right hand side vectors
     int ldab = 1 + 2*kl + ku;   //Number of rows in compressed matrix for Lapack
     int bdab = 1 + kl + ku;   //Number of rows in compressed matrix for Blas
     int ldb = n;      //size of RHS vector
     double* term1 = new double[n];
     double* term2 = new double[n];
     double* term3 = new double[n];
     double* term4 = new double[n];
     double* term5 = new double[n];



     F77Name(dgbmv)('N', n, n, klB, kuB, 1.0, B , 3, s_in, 1, 0.0, term1, 1);


     F77Name(dgbmv)('N', n, n, kl, ku, 1.0, C , bdab, v_in, 1, 0.0, term2, 1);


     F77Name(dgbmv)('N', n, n, kl, ku, 1.0, C , bdab, s_in, 1, 0.0, term3, 1);

     F77Name(dgbmv)('N', n, n, klB, kuB, 1.0, B , 3, v_in, 1, 0.0, term4, 1);
     // cout << "interior vorticity before" << endl;
     //   for (int i = 0; i < n ; i++){
     //     cout << s_in[i] << "    ";
     //   }cout << endl;
     F77Name(dgbmv)('N', n, n, kl, ku, 1.0, A_v , bdab, v_in, 1, 0.0, term5, 1);
     // cout << "interior vorticity after" << endl;
     // for (int j = 0 ; j < dNy ; j++){
     //   for (int i = 0; i < dNx ; i++){
     //     cout << term5[j + i*dNy] << "    ";
     //   }
     //   cout << "\n";
     // }

     //apply boundary conditions
     for (int j = 0 ; j < dNy ; j++){
       for (int i = 0 ; i < dNx ; i++){
         if (j == 0){         //bottom bc
             term1[j + dNy*i] -= s_bcB[i]/2.0/dy;
             term4[j + dNy*i] -= v_bcB[i]/2.0/dy;
           if (i == 0){
             term2[j + dNy*i] -= v_bcL[j]/2.0/dx;
             term3[j + dNy*i] -= s_bcL[j]/2.0/dx;
             term5[j + dNy*i] = 1.0/Re*(term5[j + dNy*i] - v_bcB[i]/dy/dy - v_bcL[j]/dx/dx);
           }
           else if(i == dNx - 1){
             term2[j + dNy*i] += v_bcR[j]/2.0/dx;
             term3[j + dNy*i] += s_bcR[j]/2.0/dx;
             term5[j + dNy*i] = 1.0/Re*(term5[j + dNy*i] - v_bcB[i]/dy/dy - v_bcR[j]/dx/dx);
           }
           else {
             term5[j + dNy*i] = 1.0/Re*(term5[j + dNy*i] - v_bcB[i]/dy/dy);
           }

         }
         else if(j == dNy - 1){      //top bc
             term1[j + dNy*i] += s_bcT[i]/2.0/dy;
             term4[j + dNy*i] += v_bcT[i]/2.0/dy;
           if (i == 0){
             term2[j + dNy*i] -= v_bcL[j]/2.0/dx;
             term3[j + dNy*i] -= s_bcL[j]/2.0/dx;
             term5[j + dNy*i] = 1.0/Re*(term5[j + dNy*i] - v_bcT[i]/dy/dy - v_bcL[j]/dx/dx);
           }
           else if(i == dNx - 1){
             term2[j + dNy*i] += v_bcR[j]/2.0/dx;
             term3[j + dNy*i] += s_bcR[j]/2.0/dx;
             term5[j + dNy*i] = 1.0/Re*(term5[j + dNy*i] - v_bcT[i]/dy/dy - v_bcR[j]/dx/dx);
           }
           else {
             term5[j + dNy*i] = 1.0/Re*(term5[j + dNy*i] - v_bcT[i]/dy/dy);
           }

         }
         else if(i == 0){    //left bc
           // term1[j + dNy*i] += dt*s_bcL[j]/2.0/dy;
           term2[j + dNy*i] -= v_bcL[j]/2.0/dx;
           term3[j + dNy*i] -= s_bcL[j]/2.0/dx;
           // term4[j + dNy*i] += dt*v_bcL[j]/2.0/dy;
           term5[j + dNy*i] = 1.0/Re*(term5[j + dNy*i] - v_bcL[j]/dx/dx);
         }
         else if(i == dNx - 1){      //right bc
           // term1[j + dNy*i] += dt*s_bcR[j]/2.0/dy;
           term2[j + dNy*i] += v_bcR[j]/2.0/dx;
           term3[j + dNy*i] += s_bcR[j]/2.0/dx;
           // term4[j + dNy*i] += dt*v_bcR[j]/2.0/dy;
           term5[j + dNy*i] = 1.0/Re*(term5[j + dNy*i] - v_bcR[j]/dx/dx);
         }
       }
     }


     for (int i = 0 ; i < n ; i++){
       term1[i] = term1[i]*term2[i];
       term3[i] = term3[i]*term4[i];
     }

     for (int i = 0 ; i < n ; i++){
       v_in[i] = v_in[i] - dt*term1[i] + dt*term3[i] - dt*term5[i];
     }

     // cout << "interior vorticity at time t + dt" << endl;
     // for (int j = 0 ; j < dNy ; j++){
     //   for (int i = 0; i < dNx ; i++){
     //     cout << v_in[j + i*dNy] << "    ";
     //   }
     //   cout << "\n";
     // }

     // delete[] term1,term2,term3,term4,term5;

 }


 void LidDrivenCavity::PoissonSolver(){

      int n = dNx*dNy;    //matrix A size
      int kl = dNy;    //Lower diagonal bandwidth
      int ku = dNy;    //Upper diagonal bandwidth
      int nrhs = 1;     //number of right hand side vectors
      int info = 1;
      int ldb = n;      //size of RHS vector
      int ldab = 1 + 2*kl + ku;    //leading size of matrix AB
      double* AB = new double[ldab*n]; //for Lapack will overwrite matrix A;
      int* ipiv = new int[n];
      double* temp = new double[n];      //for lapack will overwrite v with s;

      for (int i = 0 ; i < ldab*n; i++){
        AB[i] = A[i];
      }


      //apply boundary conditions
      for (int j = 0 ; j < dNy ; j++){
        for (int i = 0 ; i < dNx ; i++){
          if (j == 0){         //bottom bc
            if(i == 0){             //conner
              v_in[j + dNy*i] += (s_bcL[j]/dx/dx + s_bcB[i]/dy/dy);
            }
            else if(i == dNx - 1){  //conner
              v_in[j + dNy*i] += (s_bcR[j]/dx/dx + s_bcB[i]/dy/dy);
            }
            else{
              v_in[j + dNy*i] += s_bcB[i]/dy/dy;
            }
          }
          else if(j == dNy - 1){      //top bc
            if(i == 0){      //conner
              v_in[j + dNy*i] += (s_bcL[j]/dx/dx + s_bcT[i]/dy/dy);
            }
            else if(i == dNx - 1){  //conner
              v_in[j + dNy*i] += (s_bcR[j]/dx/dx + s_bcT[i]/dy/dy);
            }
            else{
              v_in[j + dNy*i] += s_bcT[i]/dy/dy;
            }
          }
          else if(i == 0){    //left bc
            v_in[j + dNy*i] += s_bcL[j]/dx/dx;
          }
          else if(i == dNx - 1){      //right bc
            v_in[j + dNy*i] += s_bcR[j]/dx/dx;
          }
        }
      }



      for (int i = 0 ; i < n ; i++){
        temp[i] = v_in[i];
      }


      // double* error_temp = new double[dNx*dNy];
      // for (int i = 0 ; i < n ; i++){
      //   error_temp[i] = s_in[i];
      // }
      // s_in_error = error_temp;
      s_in_error = s_in;

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



      F77Name(dgbsv)(n, kl, ku, nrhs, AB, ldab, ipiv, temp, n, info); //v is input and output
      // cout << "info = " << info << endl;
      // cout << "interior vorticity after" << endl;
      // for (int j = 0 ; j < dNy ; j++){
      //   for (int i = 0; i < dNx ; i++){
      //     cout << temp[j + i*dNy] << "    ";
      //   }
      //   cout << "\n";
      // }
      s_in = temp;
      // for (int j = 0 ; j < dNy ; j++){
      //   for (int i = 0; i < dNx ; i++){
      //     cout << s_in[j + i*dNy] << "    ";
      //   }
      //   cout << "\n";
      // }

      // delete[] temp;   //bug
      // delete[] AB;
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


  // double LidDrivenCavity::Error(){
  //   double er;
  //   double norm1;
  //   double norm2;
  //   for (int  i = 0 ; i < dNx*dNy ; i++){
  //     norm1 += s_in_error[i]*s_in_error[i];
  //     norm2 += s_in[i]*s_in[i];
  //   }
  //   er = fabs(sqrt(norm1) - sqrt(norm2));
  //   return er;
  // }


 double LidDrivenCavity::Error(MPI_Comm comm_cart){

   double norm1;
   double norm2;
   double er;
   double er_gather;
   for (int  i = 0 ; i < dNx*dNy ; i++){
     norm1 += s_in_error[i]*s_in_error[i];
     norm2 += s_in[i]*s_in[i];
   }
   er = fabs(sqrt(norm1) - sqrt(norm2));

   // int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
   // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count,MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)
// cout<< er <<endl;
   MPI_Reduce(&er, &er_gather, 1, MPI_DOUBLE, MPI_MAX, 0, comm_cart); //get nax error in the domain
   MPI_Bcast(&er_gather, 1, MPI_DOUBLE, 0, comm_cart); //set this value to all processors
   // cout<< er_gather <<endl;
   return er_gather;
 }







 void LidDrivenCavity::WriteToFile(){
   // int dNx = Nx - 2;
   // int dNy = Ny - 2;
   ofstream vorticityfile;
   vorticityfile.open ("vorticity.txt");
   // vorticityfile << "x   y   vorticity.\n";
   for (int j = 0 ; j < dNy ; j++){
     for (int i = 0; i < dNx ; i++){
       vorticityfile << dx*i << "   " << dy*j << "   " << v_in[j + i*dNy] << "\n";
     }
   }
   vorticityfile.close();

   ofstream streamFunctiom;
   streamFunctiom.open ("streamFunctiom.txt");
   // vorticityfile << "x   y   vorticity.\n";
   for (int j = 0 ; j < dNy ; j++){
     for (int i = 0; i < dNx ; i++){
       streamFunctiom << dx*i << "   " << dy*j << "   " << s_in[j + i*dNy] << "\n";
     }
   }
   streamFunctiom.close();
 }

 /*

 void LidDrivenCavity::test_debug(){
   int dNx = Nx-2;
   int dNy = Ny-2;
   int n = dNx*dNy;    //matrix A size
   int kl = dNy;    //Lower diagonal bandwidth
   int ku = dNy;    //Upper diagonal bandwidth
   int nrhs = 1;     //number of right hand side vectors
   int ldab = 1 + 2*kl + ku;   //Number of rows in compressed matrix for Lapack
   int bdab = 1 + kl + ku;   //Number of rows in compressed matrix for Blas
   int ldb = n;      //size of RHS vector
   cout << "for boundary vorticity" << endl;
   for (int j = 0 ; j < dNy ; j++){
     for (int i = 0; i < dNx ; i++){
       cout << v_bc[j + i*dNy] << "    ";
     }
     cout << "\n";
   }
   cout << endl;
   cout << "for interior vorticity" << endl;
   for (int j = 0 ; j < dNy ; j++){
     for (int i = 0; i < dNx ; i++){
       cout << v_in[j + i*dNy] << "    ";
     }
     cout << "\n";
   }
   cout << endl;
 }
 */


 void LidDrivenCavity::mpiSendRecive_streamf(int xs, int xd, int ys, int yd, MPI_Comm comm_cart, int cart_rank){
   double* inT_sent = new double[dNx];
   double* inB_sent = new double[dNx];
   double* inL_sent = new double[dNy];
   double* inR_sent = new double[dNy];
   for (int j = 0 ; j < dNy ; j++){
     for (int i = 0 ; i < dNx ; i++){
       if (i == 0){inL_sent[j] = s_in[j + dNy*i];}
       if (i == dNx-1){inR_sent[j] = s_in[j + dNy*i];}
       if (j == 0){inB_sent[i] = s_in[j + dNy*i];}
       if (j == dNy-1){inT_sent[i] = s_in[j + dNy*i];}
     }
   }
   // for (int i = 0 ; i < dNx ; i++){
   //   cout <<"B_r"<< s_bcB[i] << "  ";
   // }cout << endl;
   // for (int i = 0 ; i < dNx ; i++){
   //   cout <<"T_r"<< s_bcT[i] << "  ";
   // }cout << endl;
 //   if(cart_rank == 1){
 //   for (int i = 0 ; i < dNy ; i++){
 //     cout <<"cart_rank "<<cart_rank<<" L_r "<< inL_sent[i] << "  ";
 //   }cout << endl;
 // }
   // for (int i = 0 ; i < dNy ; i++){
   //   cout << xs <<"R_r"<< s_bcR[i] << "  ";
   // }cout << endl;

     MPI_Sendrecv(inR_sent, dNy, MPI_DOUBLE, xd, 1,
                     s_bcL, dNy, MPI_DOUBLE, xs, 1, comm_cart, MPI_STATUS_IGNORE);
     MPI_Sendrecv(inL_sent, dNy, MPI_DOUBLE, xs, 2,
                     s_bcR, dNy, MPI_DOUBLE, xd, 2, comm_cart, MPI_STATUS_IGNORE);
     MPI_Sendrecv(inT_sent, dNx, MPI_DOUBLE, yd, 3,
                     s_bcB, dNx, MPI_DOUBLE, ys, 3, comm_cart, MPI_STATUS_IGNORE);
     MPI_Sendrecv(inB_sent, dNx, MPI_DOUBLE, ys, 4,
                     s_bcT, dNx, MPI_DOUBLE, yd, 4, comm_cart, MPI_STATUS_IGNORE);

                  // for (int i = 0 ; i < dNx ; i++){
                  //   cout <<"B_r"<< s_bcB[i] << "  ";
                  // }cout << endl;
                  // for (int i = 0 ; i < dNx ; i++){
                  //   cout <<"T_r"<< s_bcT[i] << "  ";
                  // }cout << endl;
                  // for (int i = 0 ; i < dNy ; i++){
                  //   cout <<"L_r"<< s_bcL[i] << "  ";
                  // }cout << endl;
                //   if(cart_rank == 0){
                //     for (int i = 0 ; i < dNy ; i++){
                //       cout << xs <<"R_r"<< s_bcR[i] << "  ";
                //     }cout << endl;
                // }






                  delete[] inT_sent,inB_sent,inL_sent,inR_sent;

 }

 void LidDrivenCavity::mpiSendRecive_vorticity(int xs, int xd, int ys, int yd, MPI_Comm comm_cart){
   double* inT_sent = new double[dNx];
   double* inB_sent = new double[dNx];
   double* inL_sent = new double[dNy];
   double* inR_sent = new double[dNy];
   // for (int i = 0 ; i < dNx ; i++){
   //   cout <<"B_r"<< v_bcB[i] << "  ";
   // }cout << endl;
   // for (int i = 0 ; i < dNx ; i++){
   //   cout <<"T_r"<< v_bcT[i] << "  ";
   // }cout << endl;
   // for (int i = 0 ; i < dNy ; i++){
   //   cout <<"L_r"<< v_bcL[i] << "  ";
   // }cout << endl;
   // for (int i = 0 ; i < dNy ; i++){
   //   cout << xs <<"R_r"<< v_bcR[i] << "  ";
   // }cout << endl;
   // for (int j = 0 ; j < dNy ; j++){
   //   for (int i = 0 ; i < dNx ; i++){
   //     if (i == 0){inL_sent[j] = v_in[j + dNy*i];}
   //     if (i == dNx-1){inR_sent[j] = v_in[j + dNy*i];}
   //     if (j == 0){inB_sent[i] = v_in[j + dNy*i];}
   //     if (j == dNy-1){inT_sent[i] = v_in[j + dNy*i];}
   //   }
   // }
   // for (int i = 0 ; i < dNy ; i++){
   //   cout <<"L_s"<< inR_sent[i] << "  ";
   // }cout << endl;
     MPI_Sendrecv(inR_sent, dNy, MPI_DOUBLE, xd, 1,
                     v_bcL, dNy, MPI_DOUBLE, xs, 1, comm_cart, MPI_STATUS_IGNORE);
     MPI_Sendrecv(inL_sent, dNy, MPI_DOUBLE, xs, 2,
                     v_bcR, dNy, MPI_DOUBLE, xd, 2, comm_cart, MPI_STATUS_IGNORE);
     MPI_Sendrecv(inT_sent, dNx, MPI_DOUBLE, yd, 3,
                     v_bcB, dNx, MPI_DOUBLE, ys, 3, comm_cart, MPI_STATUS_IGNORE);
     MPI_Sendrecv(inB_sent, dNx, MPI_DOUBLE, ys, 4,
                     v_bcT, dNx, MPI_DOUBLE, yd, 4, comm_cart, MPI_STATUS_IGNORE);

                     // for (int i = 0 ; i < dNx ; i++){
                     //   cout <<"B_r"<< v_bcB[i] << "  ";
                     // }cout << endl;

                     // for (int i = 0 ; i < dNx ; i++){
                     //   cout << "from rank "<< yd <<"T_r"<< v_bcT[i] << "  ";
                     // }cout << endl;
                     // for (int i = 0 ; i < dNy ; i++){
                     //   cout <<"L_r"<< v_bcL[i] << "  ";
                     // }cout << endl;
                     // for (int i = 0 ; i < dNy ; i++){
                     //   cout << xs <<"R_r"<< v_bcR[i] << "  ";
                     // }cout << endl;
    delete[] inT_sent,inB_sent,inL_sent,inR_sent;

 }

 void LidDrivenCavity::mpiGarther(MPI_Comm comm_cart){
   // Collate local contributions
   // // double* rec_s_in = new double[dNx*dNy];
   // double* s_in_final;
   // if (comm_cart == 0){
     double* s_in_final = new double[(Nx-2)*(Ny-2)];
     // s_in_final[(Nx-2)*(Ny-2)];

   // }
   //
   // for (int coordsx=0 ; coordsx < (Nx-2) ; coordsx++){
   //   if ()
   //   MPI_Gatherv(s_in, dNy, MPI_DOUBLE, s_in_final, dNy, const int *displs, MPI_DOUBLE, 0, comm_cart)
   //
   // }
   // MPI_Gather(s_in, dNx*dNy, MPI_DOUBLE, s_in_final, dNx*dNy, MPI_DOUBLE, 0, comm_cart);
   // s_in_out = s_in_final;

   // cout << "for interior vorticity" << endl;
   // for (int j = 0 ; j < Ny-2 ; j++){
   //   for (int i = 0; i < Nx-2 ; i++){
   //     cout << s_in_out[j + i*dNy] << "    ";
   //   }
   //   cout << "\n";
   // }
   // cout << endl;


   // int MPI_Gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
   //              void *recvbuf, const int *recvcounts, const int *displs,
   //              MPI_Datatype recvtype, int root, MPI_Comm comm)

for (int j = 0 ; j < Ny-2 ; j++){
  for (int i = 0 ; i < Nx-2 ; i++){

  }
}


 }
