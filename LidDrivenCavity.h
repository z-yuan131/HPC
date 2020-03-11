#pragma once
// #include "PoissonSolver.h"
#include <string>
#include <mpi.h>


class LidDrivenCavity
{
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
    void SetSubdomainGrids(int Nx, int Ny);
    void Initialise();

    // void Initialise_top_boundary(int sizeNx, int sizeNy, double deltay);
    void CalculateVorticityBC();
    void CalculateInteriorVorticityAtTimet();
    void BuildMatrixA_B_C();
    void TimeAdvance();
    void PoissonSolver();
    double Error();
    void WriteToFile();
    void mpiSendRecive_streamf(int xr, int xd, int yr, int yd, MPI_Comm comm_cart);
    void mpiSendRecive_vorticity(int xr, int xd, int yr, int yd, MPI_Comm comm_cart);
    void test_debug();

    void Integrate();

    // friend class PoissonSolver;

    // Add any other public functions
    // There should be another function to calculate s_bc sended by another core



private:
    double* v_in = nullptr;      //array of interior vorticity
    double* v_bcL = nullptr;      //array of boundary vorticity
    double* v_bcR = nullptr;
    double* v_bcT = nullptr;
    double* v_bcB = nullptr;
    double* s_in = nullptr;      //array of interior stream function
    double* s_in_error = nullptr;
    double* s_bcT = nullptr;      //array of boundary stream function
    double* s_bcB = nullptr;
    double* s_bcL = nullptr;      //array of boundary stream function
    double* s_bcR = nullptr;
    double* A = nullptr;      //matrix(array) of linear system Ax = y
    double* A_v = nullptr;
    double* B = nullptr;      //matrix(array) for time advance
    double* C = nullptr;      //matrix(array) for time advance


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
    // friend class PoissonSolver;
};
