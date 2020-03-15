#include <iostream>
#include "LidDrivenCavity.h"
// #include "PoissonSolver.h"
#include <boost/program_options.hpp>
#include <mpi.h>

using namespace std;

int main(int argc, char **argv)
{
  //initialization
  int Nx = 30,Ny = 20;
  double Lx = 1.0,Ly = 1.0,Re = 1000.0;
  double dt = 0.0001/*Re*Lx*Ly/4/(Nx-1)/(Ny-1)*/;
  double T = 20000*dt;


  // Initialize the MPI environment
  int err = MPI_Init(&argc, &argv);
  if (err != MPI_SUCCESS) {
      cout << "Failed to initialise MPI" << endl;
      return -1;
  }

  int Px = 1,Py = 1;
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
  //      << " coords: (" << coords[0] <<","<< coords[1] << ")"<< endl;
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
  shiyu_solver->SetDomainSize(Lx,Ly);
  shiyu_solver->SetGridSize(Nx,Ny);
  shiyu_solver->Setmeshsize(Lx, Ly, Nx, Ny);
  shiyu_solver->SetPxPy(Px, Py);
  shiyu_solver->SetSubdomainGrids(Nx, Ny);
  shiyu_solver->SetFinalTime(T);
  shiyu_solver->SetTimeStep(dt);
  shiyu_solver->SetReynoldsNumber(Re);
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
  int i = 0, j;
  double error = 0.0000001;
  double er;

  do{
    j = 0;
    shiyu_solver->CalculateVorticityBC(TopBC);

    // shiyu_solver->mpiSendRecive_streamf(x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0],comm_cart);
    shiyu_solver->CalculateInteriorVorticityAtTimet();

    shiyu_solver->mpiSendRecive_vorticity(x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0],comm_cart);
    shiyu_solver->TimeAdvance();

    do{
      shiyu_solver->PoissonSolver();
      shiyu_solver->mpiSendRecive_streamf(x_rank_source[0],x_rank_dest[0],y_rank_source[0],y_rank_dest[0],comm_cart,cart_rank);
      j++;
      if (Px == 1 && Py == 1){break;}
      if (i == 0){break;}
      // cout << "j = " << j << endl;
    }while(j < 1);


    er = shiyu_solver->Error(comm_cart);
    if (cart_rank == 0){
      cout << "rank" << cart_rank << "\nTime: t = " << i*dt << "\nstep = "<< i
           << "\nnorm = "<< er << endl << endl;
    }

    i++;
  }while(dt*i < T && er > error);


    // shiyu_solver->mpiGarther(comm_cart);

    // if (cart_rank == 0){
      shiyu_solver->WriteToFile();
    // }
/**/

  // Finalize the MPI environment.
  MPI_Finalize();





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
    }*/



	return 0;
}
