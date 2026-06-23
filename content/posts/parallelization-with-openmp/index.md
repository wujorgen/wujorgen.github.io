+++ 
draft = true
date = ""
title = "Parallelization with OpenMP for Solving PDEs"
description = ""
slug = "openmp-for-pdes"
authors = []
tags = ['HPC', 'OpenMP', 'Fortran', 'PDE']
categories = ['HPC']
externalLink = ""
series = []
+++

# Background
Elliptic PDEs are notoriously expensive to solve.
Unfortunately, they commonly appear in many simulations.
This post will cover solution some iterative solution methods for elliptic PDEs.
Those are, in order:
- Jacobi Iteration
- Gauss Seidel
- Successive Over-Relaxation (SOR)
- Red-Black 

The Fortran programming language will be used instead of C or C++ due to its native support for multi-dimensional arrays.
The results will be written out to CSV files and visualized using Matplotlib in Python.

## Test Problem
The test problem will be one of the simplest examples of an elliptic PDE, Poisson's Equation:
$$
\nabla^2 \phi = f
$$

In two dimensions, this is:
$$
\frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2} = f
$$

And, it is discretized as:
$$
\frac{\phi_{i+1,j} - 2 \phi_{i, j} + \phi_{i-1, j}}{\Delta x^2} + \frac{\phi_{i,j+1} - 2 \phi_{i, j} + \phi_{i,j-1}}{\Delta x^2} =  f_{i,j}
$$

Direchlet boundary conditions will be used everywhere, as they are the simplest to implement. We will just set all points on the boundary to be zero.
The grid will be uniform in both directions so that $h=\Delta x=\Delta y$.
The source function on the right hand side is:
$$
f = 2\pi^2\sin(\pi x)\sin(\pi y)
$$

In the following, the linear system $Ax=b$ is solved, where $A$ the coefficient matrix and $a_{ij}$ are the entries of $A$. The vector $x$ is the solution to the system, and $b$ is the right hand side.

Note that for the test problem considered here, each point depends only on its immediate neighbors.
For simplicity, $\phi$ and $f$ will be stored as an $n\times n$ matrices.
In Fortran, this can be declared as:

```Fortran
real, dimension(51, 51) :: phi
real, dimension(51, 51) :: frhs
```

# Iterative Methods
Unlike _direct methods_ for solving linear systems, _iterative methods_ start with an initial guess and calculate a corrected value. This process repeats until the corrected values are sufficiently similar between iterates. It should be noted that, depending on properties of the coefficient matrix in question, iterative methods are not always guaranteed to converge.

## Jacobi Method
The Jacobi method updates all points at once based on their previous values.
$$
    x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j\neq 1} a_{ij} x_j^{(k)} \right)
$$

## Gauss Seidel
The Gauss Seidel method improves on the convergence rate of the Jacobi method by using the updated values as soon as they are available.
$$
    x_i^{(k+1)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij} x_j^{(k)} \right)
$$

## Successive Over-Relaxation (SOR)
The method of Successive Over-Relaxation aims to further improve the convergence rate of the Gauss Seidel method by "overstepping" when calculating the next value.
$$
    x_i^{(k+1)} = (1-r)x_i^{(k)} + \frac{r}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij} x_j^{(k)} \right)
$$

# Comparison of Convergence Rates
TODO

# Parallelization with Red-Black
The Red-Black scheme is more like a cross between Successive Over-Relaxation and Jacobi Method.
