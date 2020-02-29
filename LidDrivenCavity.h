#pragma once

#include <string>
using namespace std;

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
    void Initialise(int sizeNx, int sizeNy);

    // void Initialise_top_boundary(int sizeNx, int sizeNy, double deltay);
    void CalculateVorticityBC();
    void CalculateInteriorVorticityAtTimet();
    void BuildMatrixA(int Nx, int Ny);
    // void vortivityattimet();
    void Integrate();

    // Add any other public functions

private:
    double* v = nullptr;      //array of vorticity
    double* s = nullptr;      //array of stream function
    double* A = nullptr;      //matrix(array) of linear system Ax = y

    double dt;
    double T;
    int    Nx;
    int    Ny;
    double Lx;
    double Ly;
    double Re;
    double dx;
    double dy;
};
