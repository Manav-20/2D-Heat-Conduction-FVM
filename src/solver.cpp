#include "solver.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

// ---- Global parameters ----
int NX = 30;
int NY = 30;
int N  = NX * NY;

double alpha_ = 1.0;   // thermal diffusivity
double dx = 1.0, dy = 1.0;
double dt = 0.1;
int nSteps = 100;

int printInterval = 10;   // console print frequency
int writeInterval = 10;   // VTK write frequency

// ---- Functions ----
double Tbc(int i, int j) {
    if (i == 0)         return 0.0;
    if (i == NX-1)      return 100.0;
    if (j == 0)         return 0.0;
    if (j == NY-1)      return 0.0;
    return 0.0;
}

void printField(const vector<double>& T) {
    cout << fixed << setprecision(2);
    for (int j = NY-1; j >= 0; --j) {
        for (int i = 0; i < NX; ++i) {
            cout << setw(8) << T[idx(i,j)] << " ";
        }
        cout << "\n";
    }
    cout << "\n";
}

void writeVTK(const vector<double>& T, int Nx, int Ny, const string& filename) {
    ofstream vtk(filename);
    if (!vtk.is_open()) {
        cerr << "Failed to open file: " << filename << "\n";
        return;
    }
    vtk << "# vtk DataFile Version 3.0\n";
    vtk << "2D Heat Conduction\n";
    vtk << "ASCII\n";
    vtk << "DATASET STRUCTURED_POINTS\n";
    vtk << "DIMENSIONS " << Nx << " " << Ny << " 1\n";
    vtk << "ORIGIN 0 0 0\n";
    vtk << "SPACING " << dx << " " << dy << " 1\n";
    vtk << "POINT_DATA " << Nx * Ny << "\n";
    vtk << "SCALARS Temperature float 1\n";
    vtk << "LOOKUP_TABLE default\n";
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            vtk << T[idx(i,j)] << "\n";
        }
    }
    vtk.close();
    cout << "VTK written: " << filename << endl;
}

void luDecompose(const vector<vector<double>>& A,
                 vector<vector<double>>& L,
                 vector<vector<double>>& U) {
    L.assign(N, vector<double>(N, 0.0));
    U.assign(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i) L[i][i] = 1.0;

    for (int k = 0; k < N; ++k) {
        for (int j = k; j < N; ++j) {
            double sum = 0.0;
            for (int s = 0; s < k; ++s) sum += L[k][s]*U[s][j];
            U[k][j] = A[k][j] - sum;
        }
        if (fabs(U[k][k]) < 1e-14) {
            cerr << "Warning: zero pivot at k=" << k << endl;
        }
        for (int i = k+1; i < N; ++i) {
            double sum = 0.0;
            for (int s = 0; s < k; ++s) sum += L[i][s]*U[s][k];
            if (fabs(U[k][k]) < 1e-14) L[i][k] = 0.0;
            else L[i][k] = (A[i][k] - sum) / U[k][k];
        }
    }
}

vector<double> luSolve(const vector<vector<double>>& L,
                       const vector<vector<double>>& U,
                       const vector<double>& B) {
    vector<double> Y(N, 0.0), X(N, 0.0);
    for (int i = 0; i < N; ++i) {
        double sum = 0.0;
        for (int j = 0; j < i; ++j) sum += L[i][j]*Y[j];
        Y[i] = B[i] - sum;
    }
    for (int i = N-1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i+1; j < N; ++j) sum += U[i][j]*X[j];
        if (fabs(U[i][i]) < 1e-14) {
            cerr << "Warning: zero diagonal in U at i=" << i << endl;
            X[i] = 0.0;
        } else {
            X[i] = (Y[i] - sum) / U[i][i];
        }
    }
    return X;
}
