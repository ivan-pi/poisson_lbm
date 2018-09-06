# poisson_lbm

A small collection of 1D Poisson solvers based on the lattice Boltzmann method

## Description

This package contains an implementation of the one-dimensional Poisson solver described in the paper: *Chai, Z., & Shi, B. (2008). A novel lattice Boltzmann model for the Poisson equation. Applied mathematical modelling, 32(10), 2050-2058*.

## Examples

Usage examples can be found in the apps folder.

In the first example we solve the Poisson-Boltzmann equation with Debye-Huckel approximation as described in the original paper. The equation is given as <p align="center"><img src="/tex/b4a47de1c803eaf694a47143105e00f8.svg?invert_in_darkmode&sanitize=true" align=middle width=78.2611302pt height=14.202794099999998pt/></p> with the boundary conditions <p align="center"><img src="/tex/df7e036bd6299896c47c653a3ec83df1.svg?invert_in_darkmode&sanitize=true" align=middle width=149.41399934999998pt height=16.438356pt/></p>
The analytical solution of this problem is given by <p align="center"><img src="/tex/e95d7217ee0ba2eb8d99b0450be460dc.svg?invert_in_darkmode&sanitize=true" align=middle width=320.3030391pt height=40.6935375pt/></p> A plot of the analytical and numerical solutions is shown below: 

<p align="center">
  <img width="400" height="300" src="/img/example1.png">
</p>

In the second example we solve a steady-state reaction diffusion problem. For a first-order reaction in a catalyst slab we can derive the following equation:
<p align="center"><img src="/tex/8ad3d4c874e01f5947e9ad8917636ba0.svg?invert_in_darkmode&sanitize=true" align=middle width=90.61106505pt height=14.7671601pt/></p> and boundary conditions <p align="center"><img src="/tex/db89a8e5b7e96440b5d9add48c2a638d.svg?invert_in_darkmode&sanitize=true" align=middle width=171.48135014999997pt height=33.81208709999999pt/></p> 
The value Th is known as the Thiele modulus and represents the ratio between the reaction rate and the diffusion rate. An analytical solution for this problem is given by <p align="center"><img src="/tex/7b3076d4e5bae05eddbdedd6a9480aec.svg?invert_in_darkmode&sanitize=true" align=middle width=180.60250724999997pt height=38.83491479999999pt/></p> 
The symmetry (or zero-flux) boundary condition can be implemented with a second-order one-sided finite difference. A plot of the agreement between analytical and numerical solutions is shown below:
<p align="center">
  <img width="400" height="300" src="/img/example2.png">
</p>