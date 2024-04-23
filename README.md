# Comprehensive Coupled Chemistry Model (C3M)

## Description
The code is developed by the Planetary Science Laboratory at the University of Michigan, Ann Arbor. It is an atmospheric photochemistry model that can be used as a standalone photochemistry model,
and can be coupled with hydrodynamic models for chemical transport in planetary environments. It is built using the [Athena++ Magnetohydrodynamics Solver](https://github.com/PrincetonUniversity/athena), and [Cantera](https://cantera.org/), a heritage chemical kinetics code for combustion applications. C3M leverages chemical kinetics information from peer-review sources of atmospheric chemistry, open source databases, and legacy models like CalTech/JPL KINETICS.

## Installing C3M

### Pre-requisites
Following are the pre-requisites for running C3M. You can install them on your system.
```
[Cantera](https://cantera.org/)
[Athena++](https://github.com/PrincetonUniversity/athena)
[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page)
[CMake](https://cmake.org/)
```

### C3M Capabilities
```
The code currently has following capabilities:

* Solving box photochemistry, thermochemistry with ion-neutral reactions
* One dimensional gas phase neutral photochemistry and thermochemistry
* One dimensional gas phase ion-electron chemistry with electron impact processes
* Formulation of custom reaction expressions, transport variables and atmospheric opacity for photochemistry
```

### Structure of C3M
#### src
```
It contains the functions for performing radiative transfer calculations, and including custom transport, and rates for a given planet type or system of reactions. C3M leverages heritage photochemical cross sections from [MPI-Mainz](https://uv-vis-spectral-atlas-mainz.org/uvvis/cross_sections/Sulfur%20compounds/Inorganic%20S-compounds/S3.spc) database, [Leiden Observatory](https://home.strw.leidenuniv.nl/~ewine/photo/), and CalTech/JPL KINETICS.
```
#### tests
```
PhotoChemBox.cpp - It runs the box chemistry solver for a given system of chemical reactions
IonChemSolver.cpp - It solves one dimesional continuity equations for ion-neutral reactions. It includes custom reactions, and electron impact processes using second order explicit diffusion in space, and first order implicit time marching. Actinic flux calculation is based on Beer-Lambert law
1DPP.cpp - It is a modular, one dimensional photochemistry model using fully-implicit solver for a given set of chemical reactions
```
#### data
```
This directory contains important datasets corresponding to chemical reaction networks, photochemical cross sections and absorption cross sections from heritage atmospheric chemistry models.
VULCAN: This folder contains selected photochemical cross sections from Leiden Observatory stored in VULCAN exoplanet chemical kinetics package
MPIUV: This folder stores photochemical cross sections from MPI-Mainz database
planet: It contains the atmospheric structure profiles, chemical composition profiles, and electron impact rates for ionospheric chemistry (Jupiter only)
network: This folder stores Cantera and C3M compatible reaction networks in two formats (.c3m) and (.yaml)
stellar: It includes stellar intensity spectrum for various Sun and other stars
```
#### tools
```
It contains important tools to handle chemical file format, and visualization of chemical reaction networks
```
