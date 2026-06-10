+++
title = 'OpenMP for Solver Acceleration'
date = 2026-06-09T22:32:50-04:00
draft = true
summary = ""
categories = ['CFD']
tags = ['CFD', 'HPC', 'OpenMP', 'Fortran']
+++

In a previous [post]({{< relref "staggered-grid/index.md" >}}), flow over a backwards facing step was simulated using a staggered grid.
Here, that simulation is ported from Julia to Fortran.
This is done so we can use OpenMP to accelerate the solver via a red-black scheme.