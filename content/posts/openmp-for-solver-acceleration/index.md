+++ 
draft = true
date = ""
title = "OpenMP for CFD"
description = ""
slug = ""
authors = []
tags = []
categories = []
externalLink = ""
series = []
+++

Previously, flow over a backwards facing step was simulated using a [staggered grid]({{< relref "staggered-grid/index.md" >}}).
Here, that simulation is ported from Julia to Fortran.
This is done so we can use OpenMP to accelerate the solver via a red-black scheme.