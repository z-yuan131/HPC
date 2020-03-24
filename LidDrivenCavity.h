#pragma once
#include <string>
#include <mpi.h>

// class Poisson;

class LidDrivenCavity
{
// friend class Poisson;
public:
    LidDrivenCavity();    //default constructor;
    ~LidDrivenCavity();   //deconsturct default constructor

    LidDrivenCavity(double Lx, double Ly, int Nx, int Ny, int Px, int Py, double dt, double T, double Re);
    // ~LidDrivenCavity();   //deconsturct constructor


    void SetDomainSize(double xlen, double ylen);
    void SetGridSize(int nx, int ny);
    void SetTimeStep(double deltat);
    void SetFinalTime(double finalt);
    void SetReynoldsNumber(double re);
    void Setmeshsize(double lx, double ly, int nx, int ny);
    void SetPxPy(int px, int py);
    void Setparameters();
    void SetSubdomainGrids(int Nx, int Ny, int coordsx, int coordsy);
    void Initialise();

    // void Initialise_top_boundary(int sizeNx, int sizeNy, double deltay);
    void CalculateVorticityBC(int TopBC, int xs, int xd, int ys, int yd);
    void CalculateInteriorVorticityAtTimet();
    void BuildMatrixA_B_C();
    void TimeAdvance();
    // void PoissonSolver();
    double calculateprecision(MPI_Comm comm_cart);
    // double Error();
    double Error(MPI_Comm comm_cart);

    void WriteToFile(int cart_rank);
    void mpiSendRecive_streamf(int xr, int xd, int yr, int yd, MPI_Comm comm_cart,int cart_rank);
    void mpiSendRecive_vorticity(int xr, int xd, int yr, int yd, MPI_Comm comm_cart);
    void mpiGarther(MPI_Comm comm_cart, int cart_rank, int coordsx, int xs, int xd, int ys, int yd);
    void CalculateFlowVelocity();

    void test_debug();

    void Integrate();

    void printbc(int cart_rank);
    void printin(int cart_rank);


    friend class Poisson;


    // Add any other public functions
    // There should be another function to calculate s_bc sended by another core



private:
    double* v_in = nullptr;      //array of interior vorticity
    double* v_bcL = nullptr;      //array of boundary vorticity
    double* v_bcR = nullptr;
    double* v_bcT = nullptr;
    double* v_bcB = nullptr;
    double* s_in = nullptr;      //array of interior stream function
    double* s_bcT = nullptr;      //array of boundary stream function
    double* s_bcB = nullptr;
    double* s_bcL = nullptr;      //array of boundary stream function
    double* s_bcR = nullptr;
    double* A = nullptr;      //matrix(array) of linear system Ax = y
    double* A_v = nullptr;
    double* B = nullptr;      //matrix(array) for time advance
    double* C = nullptr;      //matrix(array) for time advance
    double* s_in_out = nullptr;
    double* v_in_out = nullptr;
    double* v_bc_out = nullptr;
    double* u_x = nullptr;    //flow velocity in x direction
    double* u_y = nullptr;    //flow velocity in y direction

    //calcualtion cache
    double* term1 = nullptr;
    double* term2 = nullptr;
    double* term3 = nullptr;
    double* term4 = nullptr;
    double* term5 = nullptr;
    double* term6 = nullptr;
    double* v_error = nullptr;
    double* infnorm = nullptr;
    double* inT_sent = nullptr;
    double* inB_sent = nullptr;
    double* inL_sent = nullptr;
    double* inR_sent = nullptr;






    double dt;
    double T;
    int    Nx;
    int    Ny;
    int    dNx;
    int    dNy;
    int    Px;
    int    Py;
    double Lx;
    double Ly;
    double Re;
    double dx;
    double dy;
    int    n;
    int    kl;    //Lower diagonal bandwidth
    int    ku;    //Upper diagonal bandwidth
    int    nrhs;     //number of right hand side vectors
    int    ldab;   //Number of rows in compressed matrix for Lapack
    int    bdab;   //Number of rows in compressed matrix for Blas
    int    ldb;      //size of RHS vector
    int    klB;    //Lower diagonal bandwidth
    int    kuB;    //Upper diagonal bandwidth
    double precision_gather;
    

};
