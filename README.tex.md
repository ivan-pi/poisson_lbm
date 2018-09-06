# poisson_lbm

A small collection of 1D Poisson solvers based on the lattice Boltzmann method

## Description

This package contains an implementation of the one-dimensional Poisson solver from the paper *Chai, Z., & Shi, B. (2008). A novel lattice Boltzmann model for the Poisson equation. Applied mathematical modelling, 32(10), 2050-2058*.

## Examples

Usage examples can be found in the apps folder.

In the first example we solve the Poisson-Boltzmann equation with Debye-Huckel approximation as described in the original paper. The equation is given as $$\nabla^2 u = k^2 u$$ with the boundary conditions $$u(0)=1, \quad u(1)=1.$$
The analytical solution of this problem is given by $$u(x) = \left(\frac{\exp(k) - 1}{\exp(k) - \exp(-k)}\right) \exp(-kx) + \left(\frac{1 - \exp(-k)}{\exp(k) - \exp(-k)}\right) \exp(kx).$$ A plot of the analytical and numerical solutions is shown below:


In the second example we solve a steady-state reaction diffusion problem. For a first-order reaction in a catalyst slab we can derive the following equation:
$$\nabla^2 u = \mathit{Th}^2 u$$ and boundary conditions $$u(0)=1, \quad \frac{\partial u}{\partial x}|_{x=1}=0.$$ The value \(Th\) is known as the Thiele modulus and represents the ratio between the reaction rate and the diffusion rate. An analytical solution for this problem is given by $$u(x) = \frac{cosh\left(\mathit{Th}(1-x)\right)}{\cosh(\mathit{Th})}.$$ 
The right boundary condition is implemented with a second-order one-sided finite difference. A plot of the agreement between analytical and numerical solutions is shown below: