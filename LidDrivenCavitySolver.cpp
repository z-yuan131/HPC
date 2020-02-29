#include <iostream>
using namespace std;
#include <string>
#include <boost/program_options.hpp>
// #include <mpi.h>
namespace po = boost::program_options;
#include "LidDrivenCavity.h"

#include <cstdlib>
#define F77Name(x) x##_
extern "C"{
  //LAPACK routine for solving systems of linear equations
  void F77Name(dgesv)(const int& n, const int& nrhs, const double * A,
                      const int& lda, int * ipiv, double * B,    //ipiv: vector for pivot
                      const int& ldb, int& info);
  //BLAS routine for performing matrix-vector multiplication
  void F77Name(dgemv)(const char& trans, const int& m,
                      const int& n,      const double& alpha,
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
   for (int j = 0 ; j < Ny ; j++){
     for (int i = 0 ; i < Nx ; i++){
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
     // cout << v[j+Ny*i] << " ";
   }
   // cout << "\n";
 }
}

void LidDrivenCavity::BuildMatrixA(int Nx, int Ny)
{
    int n = Nx*Ny;    //matrix A size
    int kl = Ny;    //Lower diagonal bandwidth
    int ku = Ny;    //Upper diagonal bandwidth
    int nrhs = 1;     //number of right hand side vectors
    int ldab = 1 + 2*kl + ku;   //Number of rows in compressed matrix
    int ldb = n;      //size of RHS vector
    int info;
    double* A_temp = new double[ldab*n]; //for Lapack will rewrite A matrix
    int* piv = new int[n];    //pivot data

    //populate matrix A here:
    for (int j = 0 ; j < ldab ; j++){
      for (int i = 0; i < n ; i++){
        if (j == kl+ku){        //for diagonal
          A_temp[j + i*ldab] = 4/(dx*dx + dy*dy);
        }
        else if (j == kl+ku-1 && i > 0){   //for first upper diag
          A_temp[j + i*ldab] = -1/(dx*dx);
        }
        else if (j == kl && i >= Ny){   //for second upper diag
          A_temp[j + i*ldab] = -1/(dy*dy);
        }
        else if (j == kl+ku+1 && i < n){   //for first lower diag
          A_temp[j + i*ldab] = -1/(dx*dx);
        }
        else if (j == kl+ku+kl && i < n - Ny){   //for second lower diag
          A_temp[j + i*ldab] = -1/(dy*dy);
        }
        else {
          A_temp[j + i*ldab] = 0;
        }
      cout << A_temp[j + i*ldab] << " ";}
    cout << "\n";}

}

void LidDrivenCavity::CalculateInteriorVorticityAtTimet()
{//from equation (10)

}



int main(int argc, char **argv)
{
  //initialization
  int Nx = 3,Ny = 3,Px = 1,Py = 1;
  double Lx = 1.0,Ly = 1.0,dt = 0.2,T = 1.0,Re = 100.0;
  double dx=2, dy=2;



/*

  int Nx,Ny,Px,Py;
  double Lx,Ly,dt,T ,Re;

  try {
      // Describe the options we want to have for our program
      // These are stored in a variable of type "options_description":
      po::options_description desc("Allowed options");

      // We specify the options by using the add_options() function
      // Options may take no arguments, or a single argument
      // The middle parameter specifies the type of parameter
      // The last parameter is a description displayed on the help message.
      desc.add_options()
          ("help","produce help information")
          ("Lx", po::value<double>(),"Length of the domain in the x-direction.")
          ("Ly", po::value<double>(),"Length of the domain in the y-direction.")
          ("Nx", po::value<int>(),"Number of grid points in x-direction.")
          ("Ny", po::value<int>(),"Number of grid points in y-direction.")
          ("Px", po::value<int>(),"Number of partitions in the x-direction (parallel)")
          ("Py", po::value<int>(),"Number of partitions in the y-direction (parallel)")
          ("dt", po::value<double>(),"Time step size")
          ("T", po::value<double>(),"Final time")
          ("Re",po::value<double>(),"Reynolds number")
      ;

      // These statements tell Boost program_options to parse the command-line
      // arguments given to main(), i.e. 'argc' and 'argv' in this example.
      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, desc), vm);
      po::notify(vm);

      // If the user gives the --help argument, print the help and quit.
      if (vm.count("help")) {
          cout << desc << "\n";
          return 0;
      }
int k = 0;
      // If the user entered the --compression flag, print out the value
      // given afterwards
      if (vm.count("Lx")) {
          // cout << "Value of n was set to "
          //      << vm["n"].as<int>() << ".\n";
               Lx = vm["Lx"].as<double>();k++;
      }
      else if (vm.count("Ly")) {
               Lx = vm["Ly"].as<double>();k++;
      }
      else if (vm.count("Nx")) {
               Lx = vm["Nx"].as<int>();k++;
      }
      else if (vm.count("Ny")) {
               Lx = vm["Ny"].as<int>();k++;
      }
      else if (vm.count("Px")) {
               Lx = vm["Px"].as<int>();k++;
      }
      else if (vm.count("Py")) {
               Lx = vm["Py"].as<int>();k++;
      }
      else if (vm.count("dt")) {
               Lx = vm["dt"].as<double>();k++;
      }
      else if (vm.count("T")) {
               Lx = vm["T"].as<double>();k++;
      }
      else if (vm.count("Re")) {
               Lx = vm["Re"].as<double>();k++;
      }
      else  {
          cout << "Value(s) was(were) not set, please consult --help.\n";
          return 0;
      }
  }
  // Catch any exceptions thrown (which derived from std::exception)
  // e.g. std::runtime_error, std::logic_error, etc.
  catch(exception& e) {
      cerr << "error: " << e.what() << "\n";
      return 1;
  }
  // The "..." in a catch block catches absolutely any exception thrown, no
  // matter what type it has, that hasn't been caught already by the above.
  catch(...) {
      cerr << "Exception of unknown type!\n";
  }



  */


// LidDrivenCavity LidDrivenCavity();
LidDrivenCavity*  shiyu_solver= new LidDrivenCavity(Lx,Ly,Nx,Ny,Px,Py,dt,T,Re);
shiyu_solver->SetDomainSize(Lx,Ly);
shiyu_solver->SetGridSize(Nx,Ny);
shiyu_solver->SetFinalTime(T);
shiyu_solver->SetTimeStep(dt);
shiyu_solver->SetReynoldsNumber(Re);
shiyu_solver->Setmeshsize(Lx, Ly, Nx, Ny);
shiyu_solver->Initialise(Nx, Ny);
// shiyu_solver->Initialise_top_boundary(Nx, Ny, Ly);
shiyu_solver->CalculateVorticityBC();
shiyu_solver->BuildMatrixA(Nx,Ny);
shiyu_solver->CalculateInteriorVorticityAtTimet();
// shiyu_solver->vortivityattimet();

// cout << dx<<dy<< endl;












/*
// Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    cout << "Hello world from processor " << processor_name
         << "rank " << world_rank << "out of " << world_size
         << " processors\n" << endl;

    // Finalize the MPI environment.
    MPI_Finalize();
*/




    // Create a new instance of the LidDrivenCavity class
    // LidDrivenCavity* solver = new LidDrivenCavity();

    // Configure the solver here...
    // ...
    // solver->Initialise();

    // Run the solver
    // solver->Integrate();

	return 0;
}
