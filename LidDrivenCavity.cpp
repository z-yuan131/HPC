 #include "LidDrivenCavity.h"
 using namespace std;
 #include <iostream>
 #include <string>
 #include <fstream> //write to files
 #include <algorithm>    // std::max
 #include <math.h>       /* fabs sqrt*/
 #include <boost/program_options.hpp>
 #include <mpi.h>

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

 void LidDrivenCavity::SetSubdomainGrids(int Nx, int Ny, int coordsx, int coordsy){
   if (coordsx == Px-1 && coordsy == Py-1){
     dNx  = ((Nx-2) - (Nx-2)%Px)/Px + (Nx-2)%Px;
     dNy  = ((Ny-2) - (Ny-2)%Py)/Py + (Ny-2)%Py;
   }
   else if (coordsx == Px-1){
     dNx  = ((Nx-2) - (Nx-2)%Px)/Px + (Nx-2)%Px;
     dNy  = ((Ny-2) - (Ny-2)%Py)/Py;
   }
   else if (coordsy == Py-1){
     dNx  = ((Nx-2) - (Nx-2)%Px)/Px;
     dNy  = ((Ny-2) - (Ny-2)%Py)/Py + (Ny-2)%Py;
   }
   else{
     dNx  = ((Nx-2) - (Nx-2)%Px)/Px;
     dNy  = ((Ny-2) - (Ny-2)%Py)/Py;
   }
 }

 void LidDrivenCavity::Setparameters(){
     n  = dNx*dNy;
     kl = dNy;    //Lower diagonal bandwidth
     ku = dNy;    //Upper diagonal bandwidth
     nrhs = 1;     //number of right hand side vectors
     ldab = 1 + 2*kl + ku;   //Number of rows in compressed matrix for Lapack
     bdab = 1 + kl + ku;   //Number of rows in compressed matrix for Blas
     ldb = n;      //size of RHS vector
     klB = 1;    //Lower diagonal bandwidth for B matrix
     kuB = 1;    //Upper diagonal bandwidth for B matrix

 }

 void LidDrivenCavity::Initialise()
 {
     //initialise interior streamfunction and vorticity
     double* phi = new double[dNx*dNy];
     double* omega = new double[dNx*dNy];
     for (int j = 0 ; j < dNy ; j++){
       for (int i = 0 ; i < dNx ; i++){
         phi[j + dNy*i] = 0;
         omega[j + dNy*i] = 0;
       }
     }
     //initialise boundary values of streamfunction and vorticity
     double* phi_bcL = new double[dNy];
     double* phi_bcR = new double[dNy];
     double* phi_bcB = new double[dNx];
     double* phi_bcT = new double[dNx];
     double* omega_bcL = new double[dNy];
     double* omega_bcR = new double[dNy];
     double* omega_bcB = new double[dNx];
     double* omega_bcT = new double[dNx];

     for (int i = 0 ; i < dNx ; i++){
       phi_bcT[i] = 0;           //1-2lines: top-bottom;
       phi_bcB[i] = 0;
       omega_bcT[i] = 0;           //1-2lines: top-bottom;
       omega_bcB[i] = 0;
     }
     for (int j = 0 ; j < dNy ; j++){
       phi_bcL[j] = 0;           //1-2lines: left-right;
       phi_bcR[j] = 0;
       omega_bcL[j] = 0;           //1-2lines: left-right;
       omega_bcR[j] = 0;
     }


     //update global pointers
     s_in = phi;
     v_in = omega;
     s_bcL = phi_bcL;
     s_bcR = phi_bcR;
     s_bcB = phi_bcB;
     s_bcT = phi_bcT;
     v_bcL = omega_bcL;
     v_bcR = omega_bcR;
     v_bcB = omega_bcB;
     v_bcT = omega_bcT;


     //set memory allocation for calculation and initialise them
     //improve calculation speed
     double* term1_t = new double[n];   //for equation 11: time advance
     double* term2_t = new double[n];
     double* term3_t = new double[n];
     double* term4_t = new double[n];
     double* term5_t = new double[n];
     double* term6_t = new double[n];   //for precision calculation
     double* v_error_t = new double[n]; //for error calculation
     double* infnorm_t = new double[n];
     for (int i = 0 ; i < n ; i++){
       term1_t[i] = 0;
       term2_t[i] = 0;
       term3_t[i] = 0;
       term4_t[i] = 0;
       term5_t[i] = 0;
       term6_t[i] = 0;
       v_error_t[i] = 0;
       infnorm_t[i] = 0;
     }
     term1 = term1_t;
     term2 = term2_t;
     term3 = term3_t;
     term4 = term4_t;
     term5 = term5_t;
     term6 = term6_t;
     v_error = v_error_t;
     infnorm = infnorm_t;


     double* inT_sent_t = new double[dNx];  //for sending boundaries
     double* inB_sent_t = new double[dNx];
     double* inL_sent_t = new double[dNy];
     double* inR_sent_t = new double[dNy];
     for (int i = 0 ; i < dNx ; i++){
       inT_sent_t[i] = 0;           //1-2lines: top-bottom;
       inB_sent_t[i] = 0;
     }
     for (int j = 0 ; j < dNy ; j++){
       inL_sent_t[j] = 0;           //1-2lines: left-right;
       inR_sent_t[j] = 0;
     }
     inT_sent = inT_sent_t;
     inB_sent = inB_sent_t;
     inL_sent = inL_sent_t;
     inR_sent = inR_sent_t;



     // cout << "interior streamFunction initial" << endl;
     // for (int j = 0 ; j < dNy ; j++){
     //   for (int i = 0; i < dNx ; i++){
     //     cout << s_in[j + i*dNy] << "    ";
     //   }
     //   cout << "\n";
     // }
     // cout << "interior streamFunction initial" << endl;
     // for (int j = 0 ; j < dNy ; j++){
     //   for (int i = 0; i < dNx ; i++){
     //     cout << v_in[j + i*dNy] << "    ";
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

       // for (int i = 0 ; i < dNx ; i++){
       //   cout << v_bcT[i] << "  ";
       // }cout << endl;
       // for (int i = 0 ; i < dNx ; i++){
       //     cout << v_bcB[i] << "  ";
       //   }cout << endl;
       // for (int i = 0 ; i < dNy ; i++){
       //     cout << v_bcR[i] << "  ";
       //   }cout << endl;
       // for (int i = 0 ; i < dNy ; i++){
       //     cout << v_bcL[i] << "  ";
       //   }cout << endl;

     // delete[] phi;
     // delete[] omega;
     // delete[] ohi_bc;


 }

 void LidDrivenCavity::Integrate()
 {
 }


 void LidDrivenCavity::CalculateVorticityBC(int TopBC, int xs, int xd, int ys, int yd)
 {//from equations (6)(7)(8)(9)

    for (int j = 0 ; j < dNy ; j++){
      for (int i = 0 ; i < dNx ; i++){
        if (i == 0){                  //left
          if(xs == -2){
            v_bcL[j] = (s_bcL[j] - s_in[j + dNy*i])*2.0/dx/dx;
          }
        }
        if (i == dNx-1){               //right
          if(xd == -2){
            v_bcR[j] = (s_bcR[j] - s_in[j + dNy*i])*2.0/dx/dx;
          }
        }
        if (j == 0){                  //bottom
          if(ys == -2){
            v_bcB[i] = (s_bcB[i] - s_in[j + dNy*i])*2.0/dy/dy;
          }
        }
        if (j == dNy-1){               //top
          if(yd == -2){
            v_bcT[i] = (s_bcT[i] - s_in[j + dNy*i])*2.0/dy/dy - TopBC*2.0*1.0/dy; //U = 1.0
          }
        }
      }
    }


        // for (int i = 0 ; i < dNx ; i++){
        //     cout<< v_bcT[i] << "T  ";
        //   }cout << endl;
        // for (int i = 0 ; i < dNx ; i++){
        //     cout << v_bcB[i] << "B  ";
        //   }cout << endl;
        // for (int i = 0 ; i < dNy ; i++){
        //     cout << v_bcR[i] << "R  ";
        //   }cout << endl;
        // for (int i = 0 ; i < dNy ; i++){
        //     cout << v_bcL[i] << "L  ";
        //   }cout << endl;

}


 void LidDrivenCavity::BuildMatrixA_B_C()    //dNx: points of interior domain.
 {

     // double* A_vv   = new double[bdab*n];
     double* A_vv   = new double[(1+ku)*n];

     double* A_temp = new double[ldab*n]; //for Lapack will rewrite A matrix
     double* B_temp = new double[3*n]; //B matrix uses blas routine
     double* C_temp = new double[bdab*n]; //C matrix uses blas routine

     //populate matrix A, B, C here:
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
         if (j < 1+ku){
           if (j == ku){        //for diagonal
             A_vv[j + i*(1+ku)] = 2.0/dx/dx + 2.0/dy/dy;
           }
           else if (j == ku-1 && (i%dNy)){   //for first upper diag
             A_vv[j + i*(1+ku)] = -1.0/dy/dy;
           }
           else if (j == 0 && i >= dNy){   //for second upper diag
             A_vv[j + i*(1+ku)] = -1.0/dx/dx;
           }
           // else if (j == ku+1 && ((i+1)%dNy)){   //for first lower diag
           //   A_vv[j + i*bdab] = -1.0/dy/dy;
           // }
           // else if (j == ku+kl && i < n - dNy){   //for second lower diag
           //   A_vv[j + i*bdab] = -1.0/dx/dx;
           // }
           else {
             A_vv[j + i*(1+ku)] = 0.0;
           }
         }
         //matrix C
         if (j < bdab){
           if (j == ku){        //for diagonal
             C_temp[j + i*bdab] = 0.0;
           }
           else if (j == ku-1 && (i%dNy)){   //for first upper diag
             C_temp[j + i*bdab] = 0.0;
           }
           else if (j == 0 && i >= dNy){   //for second upper diag
             C_temp[j + i*bdab] = 1.0/2.0/dx;
           }
           else if (j == ku+1 && ((i+1)%dNy)){   //for first lower diag
             C_temp[j + i*bdab] = 0.0;
           }
           else if (j == ku+kl && i < n - dNy){   //for second lower diag
             C_temp[j + i*bdab] = -1.0/2.0/dx;
           }
           else {
             C_temp[j + i*bdab] = 0.0;
           }
         }

         //matrix B
         if (j < 3){
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
     }
     A = A_temp;
     A_v = A_vv;
     // delete[] A_vv;
     // delete[] A_temp;
     for (int i = 0 ; i < n ; i++){
       for (int j = 0; j < bdab ; j++){

       }
     }

     //populate matrix B here:
     for (int i = 0 ; i < n ; i++){
       for (int j = 0; j < 3 ; j++){

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

     // F77Name(dgbmv)('N', n, n, kl, ku, 1.0, A_v, bdab, s_in, 1, 0.0, v_in, 1);
     F77Name(dsbmv)('U', n, ku, 1.0, A_v, 1+ku, s_in, 1, 0.0, v_in, 1);

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

           v_in[j + dNy*i] -= s_bcB[i]/dy/dy;
           if(i == 0){             //conner
             v_in[j + dNy*i] -= s_bcL[j]/dx/dx;
           }
           else if(i == dNx - 1){  //conner
             v_in[j + dNy*i] -= s_bcR[j]/dx/dx;
           }

         }
         else if(j == dNy - 1){      //top bc

           v_in[j + dNy*i] -= s_bcT[i]/dy/dy;
           if(i == 0){      //conner
             v_in[j + dNy*i] -= s_bcL[j]/dx/dx;
           }
           else if(i == dNx - 1){  //conner
             v_in[j + dNy*i] -= s_bcR[j]/dx/dx;
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






     F77Name(dgbmv)('N', n, n, klB, kuB, 1.0, B , 3, s_in, 1, 0.0, term1, 1);

     // cout << "interior vorticity before" << endl;
     //   for (int i = 0; i < n ; i++){
     //     cout << v_in[i] << "    ";
     //   }cout << endl;
     F77Name(dgbmv)('N', n, n, kl, ku, 1.0, C , bdab, v_in, 1, 0.0, term2, 1);
     // cout << "interior vorticity after" << endl;
     // for (int j = 0 ; j < dNy ; j++){
     //   for (int i = 0; i < dNx ; i++){
     //     cout << term2[j + i*dNy] << "    ";
     //   }
     //   cout << "\n";
     // }

     F77Name(dgbmv)('N', n, n, kl, ku, 1.0, C , bdab, s_in, 1, 0.0, term3, 1);





     F77Name(dgbmv)('N', n, n, klB, kuB, 1.0, B , 3, v_in, 1, 0.0, term4, 1);
     // cout << "interior vorticity before" << endl;
     //   for (int i = 0; i < n ; i++){
     //     cout << s_in[i] << "    ";
     //   }cout << endl;
     // F77Name(dgbmv)('N', n, n, kl, ku, 1.0, A_v , bdab, v_in, 1, 0.0, term5, 1);
     F77Name(dsbmv)('U', n, ku, 1.0, A_v, 1+ku, v_in, 1, 0.0, term5, 1);

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
             term5[j + dNy*i] -= v_bcB[i]/dy/dy;

           if (i == 0){
             term2[j + dNy*i] -= v_bcL[j]/2.0/dx;
             term3[j + dNy*i] -= s_bcL[j]/2.0/dx;
             term5[j + dNy*i] -= v_bcL[j]/dx/dx;
           }
           else if(i == dNx - 1){
             term2[j + dNy*i] += v_bcR[j]/2.0/dx;
             term3[j + dNy*i] += s_bcR[j]/2.0/dx;
             term5[j + dNy*i] -= v_bcR[j]/dx/dx;
           }

         }
         else if(j == dNy - 1){      //top bc
             term1[j + dNy*i] += s_bcT[i]/2.0/dy;
             term4[j + dNy*i] += v_bcT[i]/2.0/dy;
             term5[j + dNy*i] -= v_bcT[i]/dy/dy;

           if (i == 0){
             term2[j + dNy*i] -= v_bcL[j]/2.0/dx;
             term3[j + dNy*i] -= s_bcL[j]/2.0/dx;
             term5[j + dNy*i] -= v_bcL[j]/dx/dx;
           }
           else if(i == dNx - 1){
             term2[j + dNy*i] += v_bcR[j]/2.0/dx;
             term3[j + dNy*i] += s_bcR[j]/2.0/dx;
             term5[j + dNy*i] -= v_bcR[j]/dx/dx;
           }
           else {
           }

         }
         else if(i == 0){    //left bc
           term2[j + dNy*i] -= v_bcL[j]/2.0/dx;
           term3[j + dNy*i] -= s_bcL[j]/2.0/dx;
           term5[j + dNy*i] -= v_bcL[j]/dx/dx;
         }
         else if(i == dNx - 1){      //right bc
           term2[j + dNy*i] += v_bcR[j]/2.0/dx;
           term3[j + dNy*i] += s_bcR[j]/2.0/dx;
           term5[j + dNy*i] -= v_bcR[j]/dx/dx;
         }
       }
     }


     for (int i = 0 ; i < n ; i++){
       term1[i] = term1[i]*term2[i];
       term3[i] = term3[i]*term4[i];
     }

     for (int i = 0 ; i < n ; i++){
       v_in[i] = v_in[i] - dt*term1[i] + dt*term3[i] - dt*term5[i]/Re;
     }

     for (int i = 0 ; i < n ; i++){
       v_error[i] = - dt*term1[i] + dt*term3[i] - dt*term5[i]/Re;
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


double LidDrivenCavity::calculateprecision(MPI_Comm comm_cart){


    // F77Name(dgbmv)('N', n, n, kl, ku, 1.0, A_v , bdab, s_in, 1, 0.0, term6, 1);
    F77Name(dsbmv)('U', n, ku, 1.0, A_v, 1+ku, s_in, 1, 0.0, term6, 1);



    //apply boundary conditions
    for (int j = 0 ; j < dNy ; j++){
      for (int i = 0 ; i < dNx ; i++){
        if (j == 0){         //bottom bc
          term6[j + dNy*i] -= s_bcB[i]/dy/dy;
          if (i == 0){
            term6[j + dNy*i] -= s_bcL[j]/dx/dx;
          }
          else if(i == dNx - 1){
            term6[j + dNy*i] -= s_bcR[j]/dx/dx;
          }

        }
        else if(j == dNy - 1){      //top
          term6[j + dNy*i] -= s_bcT[i]/dy/dy;
          if (i == 0){
            term6[j + dNy*i] -= s_bcL[j]/dx/dx;
          }
          else if(i == dNx - 1){
            term6[j + dNy*i] -= s_bcR[j]/dx/dx;
          }


        }
        else if(i == 0){    //left bc
          term6[j + dNy*i] -= s_bcL[j]/dx/dx;
        }
        else if(i == dNx - 1){      //right bc
          term6[j + dNy*i] -= s_bcR[j]/dx/dx;
        }
      }
    }

    double precision;
    for (int i = 0 ; i < n ; i++){
      infnorm[i] = fabs(term6[i]-v_in[i]);
    }

    precision = infnorm[0];
    for (int i = 1 ; i < n ; i++){
      precision = fmax(precision,infnorm[i]);
    }

    MPI_Reduce(&precision, &precision_gather, 1, MPI_DOUBLE, MPI_MAX, 0, comm_cart); //get nax error in the domain
    MPI_Bcast(&precision_gather, 1, MPI_DOUBLE, 0, comm_cart); //set this value to all processors

    // delete[] term6,infnorm;
    return precision_gather;

  }


double LidDrivenCavity::Error(MPI_Comm comm_cart){

   double er;
   double er_gather;



   er = v_error[0];
   for (int i = 1 ; i < n ; i++){
     er = fmax(er,fabs(v_error[i]));
   }

   // int MPI_Bcast( void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm )
   // int MPI_Reduce(const void *sendbuf, void *recvbuf, int count,MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm)

   MPI_Reduce(&er, &er_gather, 1, MPI_DOUBLE, MPI_MAX, 0, comm_cart); //get max error in the domain
   MPI_Bcast(&er_gather, 1, MPI_DOUBLE, 0, comm_cart); //set this value to all processors

   // delete[] v_error;
   return er_gather;
 }







 void LidDrivenCavity::WriteToFile(int cart_rank){

   string rank= to_string(cart_rank);

   ofstream vorticityfile;
   vorticityfile.open ("vorticity.txt");
   // vorticityfile << "x   y   vorticity.\n";
   for (int i = 0 ; i < Nx-2 ; i++){
     for (int j = 0; j < Ny-2 ; j++){
       vorticityfile << dx*(i+1) << "   " << dy*(j+1) << "   " << v_in_out[j + i*(Ny-2)] << "\n";
     }
   }
   vorticityfile.close();

   ofstream vorticitybcfile;
   vorticitybcfile.open ("vorticity_bc.txt");
   // vorticityfile << "x   y   vorticity.\n";
   for (int j = 0 ; j < 2*(Ny-2+Nx-2) ; j++){
     vorticitybcfile << v_bc_out[j] << "\n";
   }
   vorticityfile.close();

   ofstream streamFunctiom;
   streamFunctiom.open ("streamFunctiom.txt");
   // vorticityfile << "x   y   vorticity.\n";
   for (int i = 0 ; i < Nx-2 ; i++){
     for (int j = 0; j < Ny-2 ; j++){
       streamFunctiom << dx*(i+1) << "   " << dy*(j+1) << "   " << s_in_out[j + i*(Ny-2)] << "\n";
     }
   }
   streamFunctiom.close();

   ofstream FlowVelocity;
   FlowVelocity.open ("FlowVelocity.txt");
   // vorticityfile << "x   y   vorticity.\n";
   for (int i = 0 ; i < Nx-2 ; i++){
     for (int j = 0; j < Ny-2 ; j++){
       FlowVelocity << dx*(i+1) << "   " << dy*(j+1) << "   " << u_x[j + i*(Ny-2)] << "   " << u_y[j + i*(Ny-2)] << "\n";
     }
   }
   streamFunctiom.close();
 }


 void LidDrivenCavity::mpiSendRecive_streamf(int xs, int xd, int ys, int yd, MPI_Comm comm_cart, int cart_rank){

   for (int j = 0 ; j < dNy ; j++){
     for (int i = 0 ; i < dNx ; i++){
       if (i == 0){inL_sent[j] = s_in[j + dNy*i];}
       if (i == dNx-1){inR_sent[j] = s_in[j + dNy*i];}
       if (j == 0){inB_sent[i] = s_in[j + dNy*i];}
       if (j == dNy-1){inT_sent[i] = s_in[j + dNy*i];}
     }
   }


     MPI_Sendrecv(inR_sent, dNy, MPI_DOUBLE, xd, 1,
                     s_bcL, dNy, MPI_DOUBLE, xs, 1, comm_cart, MPI_STATUS_IGNORE);
     MPI_Sendrecv(inL_sent, dNy, MPI_DOUBLE, xs, 2,
                     s_bcR, dNy, MPI_DOUBLE, xd, 2, comm_cart, MPI_STATUS_IGNORE);
     MPI_Sendrecv(inT_sent, dNx, MPI_DOUBLE, yd, 3,
                     s_bcB, dNx, MPI_DOUBLE, ys, 3, comm_cart, MPI_STATUS_IGNORE);
     MPI_Sendrecv(inB_sent, dNx, MPI_DOUBLE, ys, 4,
                     s_bcT, dNx, MPI_DOUBLE, yd, 4, comm_cart, MPI_STATUS_IGNORE);






    // delete[] inT_sent,inB_sent,inL_sent,inR_sent;

 }

 void LidDrivenCavity::mpiSendRecive_vorticity(int xs, int xd, int ys, int yd, MPI_Comm comm_cart){

   for (int j = 0 ; j < dNy ; j++){
     for (int i = 0 ; i < dNx ; i++){
       if (i == 0){inL_sent[j] = v_in[j + dNy*i];}
       if (i == dNx-1){inR_sent[j] = v_in[j + dNy*i];}
       if (j == 0){inB_sent[i] = v_in[j + dNy*i];}
       if (j == dNy-1){inT_sent[i] = v_in[j + dNy*i];}
     }
   }


   // for (int i = 0 ; i < dNx ; i++){
   //   cout <<"B_s"<< inB_sent[i] << "  ";
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
    // delete[] inT_sent,inB_sent,inL_sent,inR_sent;

 }

 void LidDrivenCavity::mpiGarther(MPI_Comm comm_cart, int cart_rank, int coordsx, int xs, int xd, int ys, int yd){
   // Collate local contributions
   // // double* rec_s_in = new double[dNx*dNy];
   // double* s_in_final;
   // if (comm_cart == 0){
     double* s_in_final = new double[(Nx-2)*(Ny-2)];
     double* v_in_final = new double[(Nx-2)*(Ny-2)];

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


   int dNy0,dNx0;

   // int MPI_Send(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
   // int MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status * status)
   for(int k = 0 ; k < dNx ; k++){
     MPI_Send(s_in+k*dNy, dNy, MPI_DOUBLE, 0, k, comm_cart);
     MPI_Send(v_in+k*dNy, dNy, MPI_DOUBLE, 0, k+Nx, comm_cart);
   }

   if (cart_rank == 0){
     int rank = 0;
     for(int i = 0 ; i < Px ; i++){
       // for(int rank = i*Py ; rank < Py+i*Py ; rank++){



         if (i == Px-1){dNx0 = dNx + (Nx-2)%Px;}
         else {dNx0 = dNx;}

         for(int k = 0 ; k < dNx0 ; k++){
           for(int rank = i*Py ; rank < Py+i*Py ; rank++){
             if (rank%Py == Py-1){dNy0 = dNy + (Ny-2)%Py;}
             else {dNy0 = dNy;}

             MPI_Recv(s_in_final+dNy*(rank%Py)+k*(Ny-2)+i*dNx*(Ny-2), dNy0, MPI_DOUBLE, rank, k, comm_cart, MPI_STATUS_IGNORE);
             MPI_Recv(v_in_final+dNy*(rank%Py)+k*(Ny-2)+i*dNx*(Ny-2), dNy0, MPI_DOUBLE, rank, k+Nx, comm_cart, MPI_STATUS_IGNORE);
           }
         }

       // }
     }
   }

   s_in_out = s_in_final;
   v_in_out = v_in_final;


 //   for (int j = 0 ; j < dNy*dNx ; j++){
 //       cout <<"rank = "<< cart_rank<< s_in[j] << "    ";
 //   }
 //   cout << endl;
 //   if (cart_rank == 0){
 //   cout << "rank = "<< cart_rank<< ",s in final" << endl;
 //   for (int j = 0 ; j < (Nx-2)*(Ny-2) ; j++){
 //       cout << s_in_out[j] << "    ";
 //   }
 //   cout << endl;
 // }

  // Gather boundary conditions to root processor
  double* v_bc_final = new double[2*(Nx-2+Ny-2)]; //store order: left-right-bottom-top
  if (xs == -2){ //left bc
    MPI_Send(v_bcL, dNy, MPI_DOUBLE, 0, 2, comm_cart);
  }

  if (cart_rank == 0){
      for(int rank = 0 ; rank < Py ; rank++){

        if (rank == Py-1){dNy0 = dNy + (Ny-2)%Py;}
        else {dNy0 = dNy;}

        MPI_Recv(v_bc_final+dNy*rank, dNy0, MPI_DOUBLE, rank, 2, comm_cart, MPI_STATUS_IGNORE);
      }
  }

  if (xd == -2){ //right bc
    MPI_Send(v_bcR, dNy, MPI_DOUBLE, 0, 3, comm_cart);
  }

  if (cart_rank == 0){
      for(int rank = 0 ; rank < Py ; rank++){

        if (rank == Py-1){dNy0 = dNy + (Ny-2)%Py;}
        else {dNy0 = dNy;}
        MPI_Recv(v_bc_final+(Ny-2)+dNy*rank, dNy0, MPI_DOUBLE, rank+Py*(Px-1), 3, comm_cart, MPI_STATUS_IGNORE);
      }
  }

  if (ys == -2){ //Bottom bc
    MPI_Send(v_bcB, dNx, MPI_DOUBLE, 0, 4, comm_cart);
  }

  if (cart_rank == 0){
      for(int rank = 0 ; rank < Py*Px-Py+1 ; rank = rank + Py){

        if (rank == Py*Px-Py){dNx0 = dNx + (Nx-2)%Px;}
        else {dNx0 = dNx;}
        MPI_Recv(v_bc_final+2*(Ny-2)+dNx*(rank/Py), dNx0, MPI_DOUBLE, rank, 4, comm_cart, MPI_STATUS_IGNORE);
      }
  }

  if (yd == -2){ //top bc
    MPI_Send(v_bcT, dNx, MPI_DOUBLE, 0, 5, comm_cart);
  }

  if (cart_rank == 0){
      for(int rank = Py-1 ; rank < Py*Px ; rank = rank + Py){

        if (rank == Py*Px-1){dNx0 = dNx + (Nx-2)%Px;}
        else {dNx0 = dNx;}
        MPI_Recv(v_bc_final+2*(Ny-2)+(Nx-2)+(rank/Py)*dNx, dNx0, MPI_DOUBLE, rank, 5, comm_cart, MPI_STATUS_IGNORE);
      }
  }
  v_bc_out = v_bc_final;

    // for (int j = 0 ; j < dNy ; j++){
    //     cout <<"rankL = "<< cart_rank<<"  "<< v_bcL[j] << "    ";
    // }
    // cout << endl;
    // for (int j = 0 ; j < dNy ; j++){
    //     cout <<"rankR = "<< cart_rank<<"  "<< v_bcR[j] << "    ";
    // }cout << endl;
    // for (int j = 0 ; j < dNx ; j++){
    //     cout <<"rankB = "<< cart_rank<<"  "<< v_bcB[j] << "    ";
    // }cout << endl;
    // for (int j = 0 ; j < dNx ; j++){
    //     cout <<"rankT = "<< cart_rank<<"  "<< v_bcT[j] << "    ";
    // }cout << endl;

  //   if (cart_rank == 0){
  //   cout << "rank = "<< cart_rank<< ",s bc final" << endl;
  //   for (int j = 0 ; j < 2*(Nx-2+Ny-2) ; j++){
  //       cout << v_bc_final[j] << "    ";
  //   }
  //   cout << endl;
  // }
 }

 void LidDrivenCavity::CalculateFlowVelocity(){
   //according to equation(3):
   //u_x = (s(j+1) - s(j-1))/2/dy;
   //u_y = (s(i+1) - s(i-1))/2/dx;
   double* u = new double[n];
   double* v = new double[n];

   F77Name(dgbmv)('N', n, n, klB, kuB, 1.0, B , 3, s_in, 1, 0.0, v, 1);
   F77Name(dgbmv)('N', n, n, kl, ku, 1.0, C , bdab, s_in, 1, 0.0, u, 1);

   u_x = u;
   u_y = v;

 }


 void LidDrivenCavity::printbc(int cart_rank){
   // cout << "interior vorticity before" << endl;
   //   for (int i = 0; i < n ; i++){
   //     cout << s_in[i] << "    ";
   //   }cout << endl;
   //
   // cout << "interior vorticity after" << endl;
   // for (int j = 0 ; j < dNy ; j++){
   //   for (int i = 0; i < dNx ; i++){
   //     cout << v_in[j + i*dNy] << "    ";
   //   }
   //   cout << "\n";
   // }
   if(cart_rank == 1){
   for (int i = 0 ; i < dNx ; i++){
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
     }cout << endl;
   }
 }

 void LidDrivenCavity::printin(int cart_rank){
   // cout << "interior vorticity before" << endl;
   //   for (int i = 0; i < n ; i++){
   //     cout << s_in[i] << "    ";
   //   }cout << endl;
   //
   if(cart_rank == 1){
   cout << "interior vorticity after" << endl;
   for (int j = 0 ; j < dNy ; j++){
     for (int i = 0; i < dNx ; i++){
       cout << v_in[j + i*dNy] << "    ";
     }
     cout << "\n";}
   }

 }
