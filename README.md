# poisson_lbm

A small collection of 1D Poisson solvers based on the lattice Boltzmann method

## Description

This package contains an implementation of the one-dimensional Poisson solver from the paper: Chai, Z., & Shi, B. (2008). A novel lattice Boltzmann model for the Poisson equation. Applied mathematical modelling, 32(10), 2050-2058.

## Examples

Usage examples can be found in the apps folder.

In the first example we solve the Poisson-Boltzmann equation with Debye-Huckel approximation as described in the original paper. The equation is given as <p align="center"><img src="/tex/b4a47de1c803eaf694a47143105e00f8.svg?invert_in_darkmode&sanitize=true" align=middle width=78.2611302pt height=14.202794099999998pt/></p> with the boundary conditions <p align="center"><img src="/tex/1b1aa94e5611fe395dbe295b8ead47a0.svg?invert_in_darkmode&sanitize=true" align=middle width=144.84777615pt height=16.438356pt/></p>. The analytical solution of this problem is given by <p align="center"><img src="/tex/700e4b01c113869a4a2bd641e5a462d7.svg?invert_in_darkmode&sanitize=true" align=middle width=410.66698859999997pt height=38.83491479999999pt/></p>. A plot of the analytical and numerical solutions is shown below:

In the second example we solve a steady-state reaction diffusion problem. For a first-order reaction in a catalyst slab we can derive the following equation:
<p align="center"><img src="/tex/8ad3d4c874e01f5947e9ad8917636ba0.svg?invert_in_darkmode&sanitize=true" align=middle width=90.61106505pt height=14.7671601pt/></p> and boundary conditions <p align="center"><img src="/tex/db89a8e5b7e96440b5d9add48c2a638d.svg?invert_in_darkmode&sanitize=true" align=middle width=171.48135014999997pt height=33.81208709999999pt/></p> The value <img src="/tex/760e1ca354518adc56df40605c196906.svg?invert_in_darkmode&sanitize=true" align=middle width=21.360426449999988pt height=22.831056599999986pt/> is known as the Thiele modulus and represents the ratio between the reaction rate and the diffusion rate. An analytical solution for this problem is given by <p align="center"><img src="/tex/17bb91f00e34af0a481d24450ad0d664.svg?invert_in_darkmode&sanitize=true" align=middle width=175.18057919999998pt height=38.83491479999999pt/></p>. The right boundary condition is implemented with a second-order one-sided finite difference. A plot of the agreement between analytical and numerical solutions is shown below: