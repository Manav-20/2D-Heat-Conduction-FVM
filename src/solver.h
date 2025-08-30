#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <string>

// Grid parameters (changeable from main if needed)
extern const int NX;
extern const int NY;
extern const int N;
extern const double dx, dy, dt, alpha_;
extern int nSteps;
extern int printInterval, writeInterval;

// Index mapping
inline int idx(int i, int j) { return i + j*NX; }

// Boundary condition
double Tbc(int i, int j);

// Printing and file writing
void printField(const std::vector<double>& T);
void writeVTK(const std::vector<double>& T, int Nx, int Ny, const std::string& filename);

// LU decomposition and solve
void luDecompose(const std::vector<std::vector<double>>& A,
                 std::vector<std::vector<double>>& L,
                 std::vector<std::vector<double>>& U);

std::vector<double> luSolve(const std::vector<std::vector<double>>& L,
                            const std::vector<std::vector<double>>& U,
                            const std::vector<double>& B);

#endif // SOLVER_H
