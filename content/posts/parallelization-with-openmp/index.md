+++ 
draft = true
date = ""
title = "Parallelization with Openmp for Solving PDEs"
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
This post will cover solution some solution methods for elliptic PDEs.
Those are, in order:
- Jacobi Iteration
- Gauss Seidel
- Successive Over-Relaxation (SOR)
- Red-Black 

The Fortran programming language will be used instead of C or C++ due to its native support for multi-dimensional arrays.

## Test Problem
The test problem will be one of the simplest examples of an elliptic PDE, Poisson's Equation:
$$
\nabla^2 \phi = f
$$

In two dimensions, it is discretized as:
$$
todo
$$

Direchlet boundary conditions will be used everywhere. Those will be implemented as follows, where $\omega$ represents the point on the boundary of the domain, and $\omega+n$ represents the $n$-th point inwards from the boundary.
$$
todo
$$

# Jacobi Iteration

# Gauss Seidel

# Successive Over-Relaxation (SOR)

# Red-Black