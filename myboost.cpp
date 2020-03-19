#include <iostream>
#include "myboost.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;
using namespace std;


void myboost::myboostfun(int argc, char **argv){

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
          exit(0);
      }
int k = 0;                        //this k is for count how many input do we have
      if (vm.count("Lx")) {       //get value from command window.
               Lx = vm["Lx"].as<double>();k++;
      }
      if (vm.count("Ly")) {
               Ly = vm["Ly"].as<double>();k++;
      }
      if (vm.count("Nx")) {
               Nx = vm["Nx"].as<int>();k++;
      }
      if (vm.count("Ny")) {
               Ny = vm["Ny"].as<int>();k++;
      }
      if (vm.count("Px")) {
               Px = vm["Px"].as<int>();k++;
      }
      if (vm.count("Py")) {
               Py = vm["Py"].as<int>();k++;
      }
      if (vm.count("dt")) {
               dt = vm["dt"].as<double>();k++;
      }
      if (vm.count("T")) {
               T = vm["T"].as<double>();k++;
      }
      if (vm.count("Re")) {
               Re = vm["Re"].as<double>();k++;
      }
      if(k != 9)  {       //if values are less than what we need, output this message.
          cout << "Value(s) was(were) not set, please consult --help.\n";
          exit(0);
          // return 0;
      }
  }
  // Catch any exceptions thrown (which derived from std::exception)
  // e.g. std::runtime_error, std::logic_error, etc.
  catch(exception& e) {
      cerr << "error: " << e.what() << "\n";
      exit(EXIT_FAILURE);
      // return 1;
  }
  // The "..." in a catch block catches absolutely any exception thrown, no
  // matter what type it has, that hasn't been caught already by the above.
  catch(...) {
      cerr << "Exception of unknown type!\n";
  }
}

int myboost::returnNx(){
  return Nx;
}
int myboost::returnNy(){
  return Ny;
}
int myboost::returnPx(){
  return Px;
}
int myboost::returnPy(){
  return Py;
}

double myboost::returnLx(){
  return Lx;
}
double myboost::returnLy(){
  return Ly;
}
double myboost::returndt(){
  return dt;
}
double myboost::returnT(){
  return T;
}
double myboost::returnRe(){
  return Re;
}
