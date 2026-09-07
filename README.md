# Computational Fluid Dynamics - Fortran Projects

## Overview
This repository contains a collection of Fortran numerical codes developed for the **Computational Fluid Dynamics (CFD)** course (Master's Degree in Aeronautical Engineering). 
The projects focus on implementing numerical methods to solve differential equations, evaluating time integration schemes, and analyzing the accuracy and stability of different computational approaches.

## Repository Structure

### 📁 Exercise 1: Numerical Integration of the Logistic Equation
* **Objective:** Solve non-linear Ordinary Differential Equations (ODEs) focusing on the Logistic Equation.
* **Methods:** Implicit Euler (1st-order accuracy) and Crank-Nicolson (2nd-order accuracy).
* **Key Features:** Error analysis, theoretical convergence rate verification, and data export for visualization.

### 📁 Exercise 2: Numerical Integration of the Linear Convection-Diffusion Equation
* **Objective:** Solve the linear convection-diffusion equation to study the physical mechanisms of advection and viscous diffusion.
* **Methods:** Time integration is performed using a three-stage third-order Runge-Kutta (RK3) scheme. Spatial discretization is handled via first-order forward (FW1), first-order backward (BW1), and second-order centered (CS2) schemes for the convective term. The diffusive term is discretized using a fourth-order centered (CS4) scheme.
* **Key Features:** Evaluates pure convection, pure diffusion, and full convection-diffusion cases using Gaussian and harmonic initial conditions. The project analyzes numerical artifacts such as the unconditional instability of FW1, the artificial dissipation of BW1, and the dispersion error of CS2. Stability limits are rigorously checked using CFL and Diffusion numbers.

### 📁 Exercise 3: 2D Incompressible Navier-Stokes Solver
* **Objective:** Solve the 2D incompressible Navier-Stokes equations within a rectangular domain.
* **Methods:** The solver utilizes a fractional-step, projection-based method. Time integration is performed using a three-stage Runge-Kutta (RK3) scheme. Spatial discretization relies on second-order finite differences on a staggered grid, while the pressure correction step uses a discrete Poisson equation solved with a preconditioned BiCGStab iterative solver.
* **Key Features:** Validates diffusive terms and physical boundaries through a Poiseuille flow simulation. Assesses numerical accuracy, energy conservation, and exact incompressibility enforcement via the Taylor-Green vortex benchmark. Explores non-linear vortex dynamics, shear layer roll-up, and subharmonic pairing through the free simulation of the Kelvin-Helmholtz instability.

---

## Prerequisites and Compilation
All codes are written in standard Fortran and can be compiled using `gfortran` (or any other Fortran compiler). 

To compile and run a specific exercise, navigate to its folder via terminal:
```bash
cd Exercise_1_Logistic_Equation
gfortran -o solver main.f90
./solver
