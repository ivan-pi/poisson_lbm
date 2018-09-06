# poisson_lbm

A small collection of 1D Poisson solvers based on the lattice Boltzmann method

## Description

This package contains an implementation of the one-dimensional Poisson solver described in the paper *Chai, Z., & Shi, B. (2008). A novel lattice Boltzmann model for the Poisson equation. Applied mathematical modelling, 32(10), 2050-2058*.

## Examples

Usage examples can be found in the apps folder.

In the first example we solve the Poisson-Boltzmann equation with Debye-Huckel approximation as described in the original paper. The equation is given as <img src="/tex/bf11dc31b1387eda6d47538750dfc8f5.svg?invert_in_darkmode&sanitize=true" align=middle width=78.26113019999998pt height=26.76175259999998pt/> with the boundary conditions <img src="/tex/c0d1b9053c9d7982696d6ef191669eaf.svg?invert_in_darkmode&sanitize=true" align=middle width=149.41399934999998pt height=24.65753399999998pt/>
The analytical solution of this problem is given by <img src="/tex/46fb971a134a0376ac37310074e91154.svg?invert_in_darkmode&sanitize=true" align=middle width=279.52033394999995pt height=37.80850590000001pt/> A plot of the analytical and numerical solutions is shown below: ![example1](/img/example1.png)


In the second example we solve a steady-state reaction diffusion problem. For a first-order reaction in a catalyst slab we can derive the following equation:
<img src="/tex/c98439b6f0ba16294ee3ce5a440f3ae0.svg?invert_in_darkmode&sanitize=true" align=middle width=90.61106505pt height=29.534320200000014pt/> and boundary conditions <img src="/tex/76602146de447e480ab83d3d9fc67ad5.svg?invert_in_darkmode&sanitize=true" align=middle width=167.9305584pt height=28.92634470000001pt/> 
The value *Th* is known as the Thiele modulus and represents the ratio between the reaction rate and the diffusion rate. An analytical solution for this problem is given by <p align="center"><img src="/tex/4f9f04585fef865177bedb4d60efe0a5.svg?invert_in_darkmode&sanitize=true" align=middle width=181.71932625pt height=38.83491479999999pt/></p> 
The symmetry (or zero-flux) boundary condition can be implemented with a second-order one-sided finite difference. A plot of the agreement between analytical and numerical solutions is shown below: ![example1](/img/example2.png)