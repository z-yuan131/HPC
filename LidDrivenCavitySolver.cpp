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

// cout << Lx <<" "<< Ly <<" "<< Nx <<" "<< Ny <<" "<< Px <<" "<< Py <<" "<< dt <<" "<< T <<" "<<Re <<" "<< endl;
// --Lx 1.0 --Ly 1.0 --Nx 10 --Ny 10 --Px 1 --Py 1 --dt 0.001 --T 1 --Re 100

// Nx = 15;Ny = 15;
// Lx = 1.0;Ly = 1.0;Re = 100.0;
// dt = 0.00001/*Re*Lx*Ly/4/(Nx-1)/(Ny-1)*/;
// T = 2*dt;
// Px = 2;Py = 1;


  // Initialize the MPI environment
  int err = MPI_Init(&argc, &argv);
  if (err != MPI_SUCCESS) {
      cout << "Failed to initialise MPI" << endl;
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
  // cout << "cart_rank: " << cart_rank
       // << " coords: (" << coords[0] <<","<< coords[1] << ")"<< endl;
  int TopBC;
  if (coords[1] == Py-1){
    TopBC = 1;
  }
  else{
    TopBC = 0;
  }


  int x_rank_source[1]; int x_rank_dest[1];
  MPI_Cart_shift(comm_cart, 0, 1, x_rank_source, x_rank_dest);
  int y_rank_source[1]; int y_rank_dest[1];
  MPI_Cart_shift(comm_cart, 1, 1, y_rank_source, y_rank_dest);
  // cout << "cart_rank = "<< cart_rank<<" x_rank_source" << x_rank_source[0]<<endl;

  // LidDrivenCavity LidDrivenCavity();
  LidDrivenCavity*  shiyu_solver= new LidDrivenCavity(Lx,Ly,Nx,Ny,Px,Py,dt,T,Re);
  Poisson* My_solver= new Poisson();
  shiyu_solver->SetDomainSize(Lx,Ly);
  shiyu_solver->SetGridSize(Nx,Ny);
  shiyu_solver->Setmeshsize(Lx, Ly, Nx, Ny);
  shiyu_solver->SetPxPy(Px, Py);
  shiyu_solver->SetSubdomainGrids(Nx, Ny, coords[0], coords[1]);
  shiyu_solver->SetFinalTime(T);
  shiyu_solver->SetTimeStep(dt);
  shiyu_solver->SetReynoldsNumber(Re);
  shiyu_solver->Setparameters();
  shiyu_solver->Initialise();
  shiyu_solver->BuildMatrixA_B_C();


  // shiyu_solver->Initialise_top_boundary(Nx, Ny, Ly);

  // PoissonSolver* Poisson= new PoissonSolver();

/*
  shiyu_solver->CalculateVorticityBC();
  shiyu_solver->CalculateInteriorVorticityAtTimet();
  shiyu_solver->mpiSendRecive_streamf(x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0],comm_cart);
  shiyu_solver->mpiSendRecive_vorticity(x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0],comm_cart);
  shiyu_solver->TimeAdvance();


  // shiyu_solver->mpiSendRecive(coords[0],coords[1],comm_cart);
  shiyu_solver->PoissonSolver();
  shiyu_solver->mpiSendRecive_streamf(x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0],comm_cart);
*/

/**/
  int i = 0,j;
  double error = 1e-4;
  double er, precision;

  My_solver->set(shiyu_solver);
  My_solver->LUfact();


  do{
    shiyu_solver->mpiSendRecive_streamf(x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0],comm_cart,cart_rank);
    shiyu_solver->CalculateVorticityBC(TopBC,x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0]);
// j = 0;

    shiyu_solver->CalculateInteriorVorticityAtTimet();

    // shiyu_solver->printin(cart_rank);
    shiyu_solver->mpiSendRecive_vorticity(x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0],comm_cart);
// shiyu_solver->printbc(cart_rank);
    shiyu_solver->TimeAdvance();


    do{
      My_solver->PoissonSolver(shiyu_solver);
      // shiyu_solver->PoissonSolver();
      if (Px == 1 && Py == 1){break;}

      shiyu_solver->mpiSendRecive_streamf(x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0],comm_cart,cart_rank);
      precision = shiyu_solver->calculateprecision(comm_cart);

      // if(cart_rank == 0){
        // cout <<"precision = " << precision<< endl;
      // }
// j++;
      if (i == 0){break;}
      // cout << "j = " << j << endl;
    }while(precision > 1e-6);




    er = shiyu_solver->Error(comm_cart);
    if (cart_rank == 0){
      cout << "rank" << cart_rank << "\nTime: t = " << i*dt << "\nstep = "<< i
           << "\nerror = "<< er << endl << endl;
    }

    i++;
  }while(dt*i < T && er > error);

    cout <<"Computation is complicated!"<<endl;
    shiyu_solver->mpiGarther(comm_cart, cart_rank, coords[0], x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0]);

    if (cart_rank == 0){
      shiyu_solver->CalculateFlowVelocity();
      shiyu_solver->WriteToFile(cart_rank);
    }

    cout << "write to file SUCCESS!" <<endl;
/**/

  // Finalize the MPI environment.
  MPI_Finalize();



	return 0;
}
