#include <iostream>
#include "LidDrivenCavity.h"
#include "Poisson.h"
#include "myboost.h"
#include <boost/program_options.hpp>
#include <mpi.h>

using namespace std;

int main(int argc, char **argv)
{
  //initialization
  int Nx,Ny,Px,Py;
  double Lx,Ly,dt,T ,Re;

  //Using boost_program_option library to get variables from command line.
  myboost* inputFun = new myboost;    //new instance input boost class "myboost"
  inputFun->myboostfun(argc,argv);    //This function get arguments from command line.
  Nx = inputFun->returnNx();          //set value of what we need.
  Ny = inputFun->returnNy();
  Px = inputFun->returnPx();
  Py = inputFun->returnPy();
  Lx = inputFun->returnLx();
  Ly = inputFun->returnLy();
  dt = inputFun->returndt();
  T  = inputFun->returnT();
  Re = inputFun->returnRe();
  inputFun->~myboost();               //eliminate instance class



  // Initialize the MPI environment
  int err = MPI_Init(&argc, &argv);
  if (err != MPI_SUCCESS) {
      cout << "Failed to initialise MPI" << endl;    //to varify -np is compatible with command line arguments --Px --Py
      return -1;
  }

  //A Cartesian grid topology is created
  int ndims = 2;                //dimensions
  int dims[ndims] = {Px, Py};   //processors in each direction
  int periods[ndims] = {0, 0};  //periodic(=1) or not(=0)
  int reorder = 1;              //allow to reorder(=1)
  MPI_Comm comm_cart;           //comm_cart is communicator with new cartesian topology (handle)
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm_cart);

  int cart_size;
  MPI_Comm_size(comm_cart, &cart_size);

  // Get the rank of the process
  int cart_rank;
  MPI_Comm_rank(comm_cart, &cart_rank);

  //the coordinates of a process
  int coords[ndims];
  MPI_Cart_coords(comm_cart, cart_rank, cart_size, coords);
  if (Px == 1 && Py == 1){
    coords[1] = 0;
  }
  // cout << "cart_rank: " << cart_rank             //to illustrate relationship of rank and coordinates
       // << " coords: (" << coords[0] <<","<< coords[1] << ")"<< endl;
  int TopBC;                                        //to identify if this rank located in top boundary of carvity
  if (coords[1] == Py-1){
    TopBC = 1;
  }
  else{
    TopBC = 0;
  }

  //get relationship between current rank and its neighboors
  int x_rank_source[1]; int x_rank_dest[1];
  MPI_Cart_shift(comm_cart, 0, 1, x_rank_source, x_rank_dest);
  int y_rank_source[1]; int y_rank_dest[1];
  MPI_Cart_shift(comm_cart, 1, 1, y_rank_source, y_rank_dest);

  //creat instance class for main solver and Poisson solver
  LidDrivenCavity*  shiyu_solver= new LidDrivenCavity(Lx,Ly,Nx,Ny,Px,Py,dt,T,Re);    //main solver
  Poisson* My_solver= new Poisson();            //Poisson Solver
  shiyu_solver->SetDomainSize(Lx,Ly);
  shiyu_solver->SetGridSize(Nx,Ny);
  shiyu_solver->Setmeshsize(Lx, Ly, Nx, Ny);
  shiyu_solver->SetPxPy(Px, Py);
  shiyu_solver->SetSubdomainGrids(Nx, Ny, coords[0], coords[1]);
  shiyu_solver->SetFinalTime(T);
  shiyu_solver->SetTimeStep(dt);
  shiyu_solver->SetReynoldsNumber(Re);
  shiyu_solver->Setparameters();              //pre-set and -allocate parameters used in the calculation
  shiyu_solver->Initialise();                 //initialise and preallocate possible vectors used in the calculation
  shiyu_solver->BuildMatrixA_B_C();           //A is for PoissonSolver, B and C is for eq. (11)


  int i = 0;
  double error = 1e-4;                        //final error control illustrate steady state.
  double er, precision;                       //precision illustrate if PoissonSolver reaches convergence

  My_solver->set(shiyu_solver);               //preset vector for Poisson Solver nad get private varibles from LidDrivenCavity class
  My_solver->LUfact();                        //pre LU factorization to parameter matrix A used in PoissonSolver.


  do{
    shiyu_solver->mpiSendRecive_streamf(x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0],comm_cart,cart_rank); //send-recive stream function.
    shiyu_solver->CalculateVorticityBC(TopBC,x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0]);                //eq. (6)-(9)

    shiyu_solver->CalculateInteriorVorticityAtTimet();                      //eq (10)

    shiyu_solver->mpiSendRecive_vorticity(x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0],comm_cart);
    shiyu_solver->TimeAdvance();                                            //eq (11)


    do{
      My_solver->PoissonSolver(shiyu_solver);                               //eq (12)
      if (Px == 1 && Py == 1){break;}                                       //if doing serial, do not need to iterate PoissonSolver

      shiyu_solver->mpiSendRecive_streamf(x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0],comm_cart,cart_rank);
      precision = shiyu_solver->calculateprecision(comm_cart);              //precision control
      // if(cart_rank == 0){        //precision can be printed on screan if needed
        // cout <<"precision = " << precision<< endl;
      // }
      if (i == 0){break;}           //for the first time loop, boundary condition is precise.
    }while(precision > 1e-6);       //precision was set to 10^-6


    er = shiyu_solver->Error(comm_cart);      //error control, decide when to stop time loop
    if (cart_rank == 0){                      //print thing on screan if needed
      cout << "rank" << cart_rank << "\nTime: t = " << i*dt << "\nstep = "<< i
           << "\nerror = "<< er << endl << endl;
    }

    i++;
  }while(dt*i < T && er > error);             //Steady T should be around 25s physical time under this error.

    shiyu_solver->CalculateFlowVelocity();    //post-process to produce velocity information

    cout <<"Computation is complicated!"<<endl;
    //gather information from other rank to rank0
    shiyu_solver->mpiGarther(comm_cart, cart_rank, x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0]);

    if (cart_rank == 0){
      shiyu_solver->WriteToFile(cart_rank);   //write information into files
    }

    cout << "write to file SUCCESS!" <<endl;


  // Finalize the MPI environment.
  MPI_Finalize();



	return 0;
}
