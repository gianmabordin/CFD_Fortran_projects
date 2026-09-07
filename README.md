# Computational Fluid Dynamics - Fortran Projects

## Overview
This repository contains a collection of Fortran numerical codes developed for the **Computational Fluid Dynamics (CFD)** course (Master's Degree in Aeronautical Engineering). 
The projects focus on implementing numerical methods to solve differential equations, evaluating time integration schemes, and analyzing the accuracy and stability of different computational approaches.

## Repository Structure

### 📁 Exercise 1: Numerical Integration of the Logistic Equation
* **Objective:** Solve non-linear Ordinary Differential Equations (ODEs) focusing on the Logistic Equation.
* **Methods:** Implicit Euler (1st-order accuracy) and Crank-Nicolson (2nd-order accuracy).
* **Key Features:** Error analysis, theoretical convergence rate verification, and data export for visualization.

---

## Prerequisites and Compilation
All codes are written in standard Fortran and can be compiled using `gfortran` (or any other Fortran compiler). 

To compile and run a specific exercise, navigate to its folder via terminal:
```bash
cd Exercise_1_Logistic_Equation
gfortran -o solver main.f90
./solver
