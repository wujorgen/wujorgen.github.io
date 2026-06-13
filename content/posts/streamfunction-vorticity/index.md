+++ 
draft = false
date = 2026-06-08T23:20:42-04:00
title = 'Streamfunction Vorticity'
description = ""
slug = ""
authors = []
tags = ['CFD', 'Julia']
categories = ['CFD']
externalLink = ""
series = ['CFD-in-Julia']
+++

The following post details a code I wrote as part of a CFD course.
It simulates flow over a backwards-facing step using the vorticity-streamfunction formulation.
Specifically, the assignment was to replicate the results from Roache and Mueller's 1970 paper "Numerical Solutions of Laminar Separated Flows".

# Vorticity-Streamfunction & Solution Procedure

The streamfunction $\psi$ is defined as:
$$
u = \frac{\partial \psi}{\partial y}, \quad v = -\frac{\partial \psi}{\partial x}
$$

The vorticity is defined as:
$$
\omega = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}
$$

The vorticity transport equation can be obtained by substituting the definition for streamfunction into the Navier-Stokes momentum equations.
$$
\frac{\partial \omega}{\partial t} + u \frac{\partial \omega}{\partial x} + v \frac{\partial \omega}{\partial y} = \frac{1}{\text{Re}}\left(\frac{\partial^2\omega}{\partial x^2}+\frac{\partial^2\omega}{\partial y^2}\right)
$$

The solution procedure for the vorticity-streamfunction formulation (VSF) of the Navier-Stokes Equations (NSE) is as follows.

- Calculate $\omega$ from a velocity field $u$, $v$.
- March $\omega$ forward in time using the transport equation.
    - A first-order explicit Euler timestep ($\omega^{next}=\omega + \Delta t (\partial\omega/\partial t)$ and centered finite differences will be used for simplicity. Boundary conditions for $\omega$ will be enforced during this step.
- Solve for the streamfunction $\psi$ from the transported vorticity.
    - This requires the solution of an elliptic PDE ($\partial_{xx}\psi + \partial_{yy}\psi = -\omega$). The method of successive over-relaxation will be used. Boundary conditions for $\psi$ will be enforced during this step.
- Recover the velocity field $u$, $v$ from the streamfunction $\psi$.
    - Velocity boundary conditions are enforced during this step, and before the calculation of the vorticity in Step 1.
- The previous steps are repeated until a steady-state solution is reached. Afterwards, the pressure can be recovered (by solution of a Poisson equation).


# Boundary Conditions
In their work, Roache and Mueller simulated flow over the backwards-facing step with a free surface at the top of the domain.
The stated reason is that they are looking for the least restrictive and most natural boundary conditions.
The exact quote is: *"The boundary condition at the upper boundary or lid was somewhat disappointing. It had been desired to model the backstep with no upper boundary (free-flight case), but numerical experimentation did not disclose a method of allowing inflow through the mesh at the lid without destabilizing the solution."*
They also do not clearly specify the inlet boundary condition beyond stating that it is arrived at by integrating a Pohlhausen equation.
For simplicity (and maybe greater disappointment), this replication will alter the boundary conditions slightly to better constrain the problem.

The inlet is assumed to be a fully developed Poiseuille flow with a parabolic velocity profile. (As the simulation is 2D, it is assumed that the velocity profile corresponds to a fully developed flow between parallel plates.)
This allows for a uniform pressure distribution to be assumed at the inlet, and for this inlet pressure to be treated as a Dirichlet boundary condition.
This also has the benefit of simplifying the pressure recovery, which will be discussed later.
Finally, the top of the domain is treated with a no-slip condition rather than attempting to simulate a free surface.

It must be noted that these alterations will result in differences between Roache and Mueller's work and the results herein - the plots will not match exactly.
The location and size of recirculation zones behind the step will differ, as well as the exact contours for vorticity, pressure, and the streamfunction.
Direct comparisons will not be possible; instead, qualitative comparisons will be done to verify the code.

The following subsections discuss each of the boundary conditions necessary to define the simulation.

## Velocity
As stated previously, the inlet is placed on the left side of the domain, and a fully developed Poiseuille flow is assumed.
The velocity profile for a Poiseuille flow between two parallel plates, with the plates located at $y=h$ and $y=-h$ (so that the centerline of the flow path is at $y=0$), can be obtained with the following equation.
In the simulation, the inlet will not be between $\pm h$, but it is a relatively simple fix to translate the profile vertically.
<!--\begin{equation} \label{eqn:bc_u_inlet}-->
<div>
$$
u(y)\big|_{inlet} = u_{max}\bigg( 1 - \big(\frac{y}{h}\big)^2 \bigg)
$$
</div>
<!--\end{equation}-->
For any fully developed Poiseuille flow, the maximum velocity occurs at the centerline of the flow path.
There exist relations between the average and maximum velocity of the velocity profile, but they are not necessary here as the inlet boundary condition will be specified by $u_{max}$.

The outlet is on the right side of the domain.
Unlike the inlet, no velocity profile is assumed.
The only condition to be enforced at the outlet is that there is no velocity gradient, so that:
<!--\begin{equation} \label{eqn:bc_u_out}-->
<div>
$$
\frac{\partial u}{\partial x}\bigg|_{outlet}=0, \quad \frac{\partial v}{\partial x}\bigg|_{outlet}=0
$$
</div>
<!--\end{equation}-->

Impermeable walls with no-slip boundary conditions are assumed along the top and bottom of the domain.
<!--\begin{equation} \label{eqn:bc_u_topbot}-->
<div>
$$
u(x)\big|_{top}=0, \quad u(x)\big|_{bottom}=0
$$
</div>
<!--\end{equation}-->

## Streamfunction
The boundary conditions for the streamfunction must be calculated from the velocity boundary conditions.
As discussed previously, there are impermeability and no-slip boundaries at the top and bottom of the domain.
These are written as follows (including the definition of the streamfunction):
$$
v=-\frac{\psi}{\partial x} = 0, \quad u = \frac{\psi}{\partial y} = 0,
$$
For the top and bottom of the domain, it can be found that $\psi=C$ where $C$ is an integration constant.

To determine the value of this integration constant, the inlet to the domain must be considered.
The inlet velocity profile $u(y)$ is known, which means that $\partial \psi / \partial y$ is also known.
The streamfunction inlet boundary condition can then be obtained by integrating the $u$-velocity profile in the $y$ direction.
To set the integration constant, it will be assumed that the streamfunction is zero along the top of the domain, so that $\psi\big|_{top}=0$.
Note that the same $\pm h$ is assumed as in Equation \ref{eqn:bc_u_inlet}, as it is trivial to shift the profile up and down in code.
<!--\begin{equation} \label{eqn:bc_stream_inlet}-->
<span style="display: block; overflow-x: auto; padding-bottom: 10px;">
$$
\psi(y)\big|_{inlet} = u_{max} \bigg( y - \frac{y^3}{3\Delta y^2} \bigg)+C, \quad C = -\frac{2u_{max}\Delta y}{3}
$$
</span>
<!--\end{equation}-->
With this result in mind, the streamfunction is set to $0$ along the top of the domain and $-\frac{2u_{max}\Delta y}{3}$ along the bottom of the domain.

For the outflow condition, the streamfunction is simply extrapolated.
<!--\begin{equation} \label{eqn:bc_stream_outlet}-->
<div>
$$
\psi\big|_{N_x} = 2 \psi\big|_{N_x-1} - \psi\big|_{N_x-2}
$$
</div>
<!--\end{equation}-->

## Vorticity
As with the velocity and streamfunction, the vorticity at the inlet is known due to the assumption of a fully developed Poiseuille flow.
With the assumption that $\partial v / \partial x=0$ at the inlet, it is simply:
<!--\begin{equation} \label{eqn:bc_w_inlet}-->
<div>
$$
\omega(y)\big|_{inlet} = -\frac{\partial u}{\partial y} = \frac{2u_{max}y}{\Delta y^2}
$$
</div>
<!--\end{equation}-->

Unlike the streamfunction, the vorticity does not need to be held at a constant value along the top and bottom of the domain.
Relations such as Thom's method or Jensen's method must be used to determine the vorticity at the wall from the streamfunction.
Roache and Mueller use Thom's method, which is presented below.
Note that the $w$ subscript denotes a point on the wall, with $w+1$ indicating an interior point.
<!--\begin{equation} \label{eqn:bc_w_wall}-->
<div>
$$
\omega_w = -2(\psi_{w+1}-\psi_w)/\Delta n^2
$$
</div>
<!--\end{equation}-->

For the outlet, a simple no-gradient assumption is made.
<!--\begin{equation} \label{eqn:bc_w_out}-->
<div>
$$
\frac{\partial \omega}{\partial x} \bigg|_{outlet} = 0
$$
</div>
<!--\end{equation}-->

# Pressure Recovery
To recover the pressure after a steady-state solution has been obtained, the divergence of the momentum equations is taken, and the resulting Poisson equation is solved for pressure.
The pressure terms in the momentum equations are isolated as follows.
Note that, since the pressure is only recovered after a steady-state solution has been reached, $\partial u/\partial t$ and $\partial v/\partial t$ are $0$.
<!--\begin{equation} \label{eqn:dpdn}-->
$$
\begin{split}
    \frac{\partial p}{\partial x} = \frac{1}{Re}\bigg(\frac{\partial^2u}{\partial x^2}+\frac{\partial^2u}{\partial y^2}\bigg) - u\frac{\partial u}{\partial x} - v\frac{\partial u}{\partial y} \equiv f \\\\
    \frac{\partial p}{\partial y} = \frac{1}{Re}\bigg(\frac{\partial^2v}{\partial x^2}+\frac{\partial^2v}{\partial y^2}\bigg) - u\frac{\partial v}{\partial x} - v\frac{\partial v}{\partial y} \equiv g \\\\
\end{split}
$$
<!--\end{equation}-->
The divergence of the isolated momentum terms results in a Poisson equation for the pressure.
$$
\nabla^2 p = \frac{\partial}{\partial x}f + \frac{\partial}{\partial y}g
$$
In the following, $\partial_x\equiv \frac{\partial}{\partial x}$, and similarly for $y$.
This is done to enhance readability (and for ease of typing).
$$
\partial_x f = \frac{1}{Re}(\partial_{xxx}u + \partial_{yyx}u) - \partial_x(u\partial_x u + v\partial_y u)
$$
$$
\partial_x g = \frac{1}{Re}(\partial_{xxy}v + \partial_{yyy}v) - \partial_y(u\partial_x v + v\partial_y v)
$$
Plugging these back into the pressure Poisson equation (note the terms that equal zero due to continuity):
<span style="display: block; overflow-x: auto; padding-bottom: 10px;">
$$
\nabla^2 p = \frac{1}{Re}\big( \cancel{\partial_{xxx}u+\partial_{xxy}v} + \cancel{\partial_{yyx}u+\partial_{yyy}v} \big) - \partial_x(u\partial_x u + v\partial_y u)- \partial_y(u\partial_x v + v\partial_y v)
$$
$$
\nabla^2 p = -\partial_xu\partial_xu - u\partial_{xx}u - \partial_xv\partial_yu - v\partial_{xy}u - \partial_yu\partial_xv - u\partial_{xy}v - \partial_yv\partial_yv - v\partial_{yy}v
$$
$$
\nabla^2 p = -(\partial_x u)^2 - 2\partial_xv\partial_yu - (\partial_yv)^2 -u\partial_x\cancel{(\partial_xu+\partial_yv)} - v\partial_y\cancel{(\partial_xu+\partial_yv)}
$$
<!--\begin{equation} \label{eqn:pressure_poisson}-->
$$
\nabla^2 p = -(\partial_x u)^2 - 2\partial_xv\partial_yu - (\partial_yv)^2
$$
</span>
<!--\end{equation}-->

The equations for $\partial p / \partial x$ and $\partial p / \partial y$ can be used to set Neumann boundary conditions for the equation shown above, specifically along the top and bottom of the domain.
The outlet of the domain is assumed to have a simple Neumann boundary condition so that $\partial p/\partial x = 0$.
Note that if Neumann boundary conditions are assumed everywhere, this system will be indeterminate and will not converge, as the solution is unique up to an additive constant.
As described previously, the inlet pressure boundary was set to a Dirichlet boundary condition to simplify pressure recovery by eliminating this issue.
The pressure at the inlet is then set to zero.

# Simulation Results
The simulation domain is eight units wide and one unit tall.
Each unit is discretized with 61 grid points.
The backwards-facing step was placed in the bottom-left corner of the domain and has a height of half a unit and a length of one unit.
The maximum velocity of the inlet is set to $1.0$, and, as the non-dimensional form of the Navier-Stokes equations is used, the Reynolds number is directly set.
As transient equations were described previously, a steady-state solution was obtained by marching in time.
Convergence was determined when the $L2$-norm of the difference between iterations for both the vorticity and streamfunction decreased below a threshold of 1E-6.

Below are streamfunction contour plots for a Reynolds number of 100, 250, and 400.

![](streamfunction100.png)
![](streamfunction250.png)
![](streamfunction400.png)

Below are vorticity plots for a Reynolds number of 100, 250, and 400.

![](vorticity100.png)
![](vorticity250.png)
![](vorticity400.png)

And last but not least, below are pressure contour plots for Reynolds numbers of 100, 250, and 400.

![](pressure100.png)
![](pressure250.png)
![](pressure400.png)

## Effect of Reynolds Number
<!--A contour plot of the streamfunction for multiple Reynolds numbers is shown in Figure \ref{fig:streamfunctions}.-->
<!--Similarly, a contour plot of the vorticity and pressure fields are shown in Figure \ref{fig:vorticity} and \ref{fig:pressure}.-->
The following effects can be observed as the Reynolds number increases:
- The recirculation zone behind the step increases in size.
- The point at which the flow reattaches to the lower edge of the domain moves further downstream. This is most easily observed in the plots of the streamfunction and the pressure field.
- The region of negative vorticity at the sharp corner of the step extends further to the right and down. This seems to be related to Roache and Mueller's observation of how the flow separation point moves down the vertical face of the step as the Reynolds number increases.
- The region of negative pressure at the sharp corner of the step increases in size but decreases in magnitude. This may be a result of decreased viscous stresses, and it is possible that the pressure spike will increase in magnitude again for higher Reynolds numbers.

# Source Code

The code can be viewed at this [link](https://gist.github.com/wujorgen/ddc48c89a7430f36e33262e4237c5e87).