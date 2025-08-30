#include "solver.h"
#include <iostream>
#include <vector>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[]) {

  // Allow user to override defaults via command line
    if (argc >= 3) { NX = atoi(argv[1]); NY = atoi(argv[2]); }
    if (argc >= 4) alpha_ = atof(argv[3]);
    if (argc >= 5) dx = dy = atof(argv[4]);
    if (argc >= 6) dt = atof(argv[5]);
    if (argc >= 7) nSteps = atoi(argv[6]);

    N = NX * NY; // update dependent variable

    cout << "Simulation parameters:\n";
    cout << "Grid: " << NX << " x " << NY << "\n";
    cout << "alpha = " << alpha_ << ", dx=dy=" << dx
         << ", dt=" << dt << ", steps=" << nSteps << "\n\n";
    // Precompute directional Fo
    const double Fox = alpha_ * dt / (dx*dx);
    const double Foy = alpha_ * dt / (dy*dy);

    // Build system matrix A
    vector<vector<double>> A(N, vector<double>(N, 0.0));
    for (int j = 0; j < NY; ++j) {
        for (int i = 0; i < NX; ++i) {
            int p = idx(i,j);
            bool isBnd = (i==0 || i==NX-1 || j==0 || j==NY-1);
            if (isBnd) {
                A[p][p] = 1.0;
            } else {
                double aW = Fox, aE = Fox, aS = Foy, aN = Foy;
                double aP = 1.0 + 2.0*(Fox + Foy);
                A[p][p]               = aP;
                A[p][idx(i-1,j)]      = -aW;
                A[p][idx(i+1,j)]      = -aE;
                A[p][idx(i,  j-1)]    = -aS;
                A[p][idx(i,  j+1)]    = -aN;
            }
        }
    }

    // LU decomposition
    vector<vector<double>> L, U;
    luDecompose(A, L, U);

    // Initial condition
    vector<double> T(N, 0.0);
    for (int j=0;j<NY;++j)
        for (int i=0;i<NX;++i)
            if (i==0||i==NX-1||j==0||j==NY-1)
                T[idx(i,j)] = Tbc(i,j);

    cout << "Initial temperature field (t=0):\n";
    printField(T);

    // Time marching
    for (int n = 0; n < nSteps; ++n) {
        vector<double> B(N, 0.0);
        for (int j=0;j<NY;++j) {
            for (int i=0;i<NX;++i) {
                int p = idx(i,j);
                if (i==0||i==NX-1||j==0||j==NY-1)
                    B[p] = Tbc(i,j);
                else
                    B[p] = T[p];
            }
        }

        vector<double> Tnew = luSolve(L, U, B);

        for (int j=0;j<NY;++j)
            for (int i=0;i<NX;++i)
                if (i==0||i==NX-1||j==0||j==NY-1)
                    Tnew[idx(i,j)] = Tbc(i,j);

        T.swap(Tnew);

        if (printInterval > 0 && ((n+1)%printInterval==0)) {
            cout << "Step " << n+1 << ", t=" << (n+1)*dt << "\n";
            printField(T);
        }
        if (writeInterval > 0 && ((n+1)%writeInterval==0)) {
            string fname = "../results/case_" + to_string(n+1) + ".vtk";
            writeVTK(T, NX, NY, fname);
        }
    }

    // final write
    writeVTK(T, NX, NY, "../results/case1.vtk");
    cout << "Simulation finished.\n";
    return 0;
}
