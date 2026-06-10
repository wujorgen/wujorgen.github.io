+++
title = 'Staggered Grid'
date = 2026-06-09T22:25:25-04:00
draft = true
summary = ""
categories = ['CFD']
tags = ['CFD', 'Julia']
+++

In a previous [post]({{< relref "streamfunction-vorticity/index.md" >}}), the vorticity-streamfunction formulation of the Navier-Stokes Equations was used to simulate flow over a backwards facing step.
This post will simulate flow over the backwards facing step using the primitive variable approach and a staggered grid.