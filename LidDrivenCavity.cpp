 #include "LidDrivenCavity.h"

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
   dx = lx/nx;
   dy = ly/ny;
 }

 void LidDrivenCavity::Initialise(int sizeNx, int sizeNy)
 {
     double* phi = new double[sizeNx*sizeNy];
     double* omega = new double[sizeNx*sizeNy];
     for (int i = 0 ; i < sizeNx ; i++){
       for (int j = 0 ; j < sizeNy ; j++){
         phi[j + sizeNy*i] = 0;
         omega[j + sizeNy*i] = 0;
       }
     }
     // cout << phi[2] << endl;
     s = phi;
     v = omega;
     delete[] phi;
     delete[] omega;
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


 void LidDrivenCavity::CalculateVorticityBC()
 {//from equations (6)(7)(8)(9)
    for (int i = 0 ; i < Nx ; i++){
      for (int j = 0 ; j < Ny ; j++){
        if (i == 0){                  //left
            v[j + Ny*i] = (s[j + Ny*0] - s[j + Ny*1])*2/dx/dx;
        }
        if (i == Nx-1){               //right
            v[j + Ny*i] = (s[j + Ny*(Nx-1)] - s[j + Ny*(Nx-2)])*2/dx/dx;
        }
        if (j == 0){                  //bottom
            v[j + Ny*i] = (s[0 + Ny*i] - s[1 + Ny*i])*2/dy/dy;
        }
        if (j == Ny-1){               //top
            v[j + Ny*i] = (s[Ny-1 + Ny*i] - s[Ny-2 + Ny*i])*2/dy/dy - 2*1.0/dy; //U = 1.0
        }
      }cout << v[j+Ny*i] << " ";
    }cout << "/n";
 }

 void LidDrivenCavity::CalculateInteriorVorticityAtt()
 {

 }
