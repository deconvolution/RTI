# RTI

[![Runtests](https://github.com/deconvolution/RTI/actions/workflows/Runtests.yml/badge.svg?branch=main)](https://github.com/deconvolution/RTI/actions/workflows/Runtests.yml)
[![Documentation](https://github.com/deconvolution/RTI/actions/workflows/Documenter.yml/badge.svg?branch=main)](https://github.com/deconvolution/RTI/actions/workflows/Documenter.yml)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://deconvolution.github.io/RTI/dev/)

A Julia package for Ray Theory Inversion (RTI).
This package contains the forward solver and adjoint solver for the eikonal eqution. With both simulations, the update direction of the velocity (gradient) can be obtained. The optimization is done by l-BFGS.
## Contents
* [Target](#Target)
* [Highlights](#Highlights)
* [Installation](#Installation)
* [Dependencies](#Dependencies)
* [Usage](#Usage)
## Target
The target of this package is to create an easy-to-use package for travel time inversion.
## Highlights
Some of the highlights are:
- Shots are parallelized.
- The solvers are based on the fast-sweeping method.
- No Frechet derivatives are explicitly constructed. It is suitable for big data.
## Installation
In an IDE with Julia installed (e.g., [Atom](https://atom.io/)), one can install this package by
```julia
julia> ]
(@v1.6) pkg> add https://github.com/deconvolution/RTI
```
You can test whether it works on your system with
```julia
julia> ]
(@v1.6) pkg> test RTI
```
and use it with
```julia
julia> using RTI
```
## Other packages needed
- [MATLAB](https://github.com/JuliaInterop/MATLAB.jl)
- [Optim](https://github.com/JuliaNLSolvers/Optim.jl)
- [LineSearches](https://github.com/JuliaNLSolvers/LineSearches.jl)

## Usage
See [documentation](https://deconvolution.github.io/RTI/dev/).
