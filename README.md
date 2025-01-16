# Equilibrium Model With Anisotropy

This repository contains code for modeling the physics of magnetic
nanoparticles for the usage in model-based magnetic particle imaging (MPI).
In particular, the equilibrium model with anisotropy is in the focus of this code repository.

The methods corresponding to this code are described in the associated publication

M. Maass, T. Kluth, C. Droigk, H. Albers, K. Scheffler, A. Mertins, and T. Knopp. Equilibrium Model with Anisotropy for Model-Based Reconstruction in Magnetic Particle Imaging, IEEE Transactions on Computational Imaging, vol. 10, pp. 1588-1601, 2024.

Update from 09/01/2025:
Please note that a bugfix has been made in systemMatrix.jl.  The version corresponding to the article is now labeled “Version_Article” in Git.

## Installation

In order to use this code one first has to download [Julia](https://julialang.org/) (version 1.10 or later), clone this repository and navigate to the folder in the command line. The example scripts automatically activate the environment and install all necessary packages. This will take several minutes before the actual code is run, since all packages are precompiled during installation.

## Execution
We recommend to start Julia with multiple threads, since some computations otherwise last to long. For instance, if your CPU has 8 cores, you can start Julia with
```
julia -t 8
```
Then you can start the central script of this repository
```julia
include("simulation.jl")
```
which will calculate MPI system matrices for various models and compare them to a measured system matrix. Afterwards the system matrices an be used to perform reconstruction on various particle
phantoms by invoking
```julia
include("reconstruction.jl")
```
Both scripts will generate images in the folder `img`. 

In addition to these two scripts, there are the following scripts:
* `accuracy1D.jl`: Investigates the accuracy of the different models based on a simple 1D excitation.
* `accuracy2D.jl`: Investigates the accuracy of the different models based on a 2D Lissajous excitation.
* `order.jl`: Investigates the series truncation index within the equilibrium model with anisotropy.
These scripts have a very long computation time and should only be run in order to reproduce the figures from the associated publication.

## Open MPI Data

The measurement data associated to this project is about 80 MB large and will be downloaded and stored automatically, when the code is executed for the first time.
It is published under a [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/legalcode) license and can be found here:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10646064.svg)](https://doi.org/10.5281/zenodo.10646064)
