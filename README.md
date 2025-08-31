# 2D Heat Conduction using Finite Volume Method (FVM)

This project implements a 2D transient heat conduction solver using the Finite Volume Method (FVM) in C++.  
It demonstrates how numerical methods can be used to solve partial differential equations (PDEs) and visualize the results in ParaView via VTK output files.  

---

## Problem Description

We solve the 2D transient heat conduction equation:

∂T/∂t = α ( ∂²T/∂x² + ∂²T/∂y² )


Where:
- `T(x, y, t)` = temperature field  
- `α` = thermal diffusivity  
- The domain is discretized into a uniform structured grid with spacing `Δx`, `Δy`.  

Time marching is performed with an implicit scheme, which is unconditionally stable.

---

## Features
- Finite Volume Method (FVM) based discretization  
- Implicit time integration (better stability)  
- Flexible grid size (NX × NY) and time step  
- User can set simulation parameters (grid size, diffusivity, Δx, Δt, steps)  
- Exports results to VTK format for visualization in ParaView  
- Clean modular structure (`main.cpp`, `solver.cpp`, `solver.h`)  

---

## Project Structure

```
2D-Heat-Conduction-FVM/
│── src/ # C++ source files
│ ├── main.cpp # Entry point, user input, calls solver
│ ├── solver.cpp # Numerical solver implementation
│ └── solver.h # Header file for solver
│── results/ # Example outputs
│ ├── case1.vtk # VTK file for visualization in ParaView
│── README.md # Project documentation
│── Makefile # (optional) for easy build
```


---

## Simulation Parameters

Default values (can be overridden by user):

| Parameter | Default | Description |
|-----------|---------|-------------|
| `NX` | 30 | Number of grid cells in X |
| `NY` | 30 | Number of grid cells in Y |
| `alpha` | 1.0 | Thermal diffusivity |
| `dx, dy` | 1.0 | Grid spacing |
| `dt` | 0.1 | Time step size |
| `nSteps` | 100 | Number of iterations |

---

## Usage

### 1. Compile
```bash
cd src
g++ main.cpp solver.cpp -o heat2d
```
### 2. Run with defaults
``` bash
./heat2d
```
### 3. Run with custom parameters
```bash
./heat2d NX NY alpha dx dt nSteps
```
Example:

```bash
./heat2d 50 50 0.5 0.5 0.05 200
```
### 4. Output
The solver writes results in VTK format inside results/:

```bash
results/case1.vtk
