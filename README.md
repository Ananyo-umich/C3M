# Comprehensive Coupled Chemistry Model (C3M)

## Description
The code is developed by the Planetary Science Laboratory at the University of Michigan, Ann Arbor. It is an atmospheric photochemistry model that can be used as a standalone photochemistry model,
and can be coupled with hydrodynamic models for chemical transport in planetary environments. It is built using the Athena++ Magnetohydrodynamics Solver (https://github.com/PrincetonUniversity/athena), and Cantera (https://cantera.org/), a heritage chemical kinetics code for combustion applications. C3M leverages chemical kinetics information from peer-review sources of atmospheric chemistry, open source databases, and legacy models like CalTech/JPL KINETICS.

Following are the pre-requisites for running C3M. You can install them on your system.
'''
Cantera
Athena++
Eigen
CMake
'''

The code currently has following capabilities:

i) Solving zero dimensional photochemistry, thermochemistry with electron impact processes
ii) One dimensional gas phase neutral photochemistry and thermochemistry
iii) One dimensional gas phase ion-electron chemistry with electron impact processes
iv) Formulation of custom reaction expressions, transport variables and atmospheric opacity for photochemistry


With applications to Earth, Venus, Jupiter ionosphere



Code Structure
-> src
   -- It contains the functions for performing atmospheric chemistry and enhancing code modularity
-> tests
   -- Contains file for running test case, and model validation
   -- BoxChemSolver.cpp - Code for zero dimensional chemistry
   -- ChemSolver.cpp - One dimensional photochemistry model
   -- IonChemSolver.cpp - One dimensional photochemistry model with provision for ion impact process in atmosphere introduced artificially
-> data
   -- This directory contains important datasets corresponding to chemical reaction networks, photochemical cross sections and absorption cross sections from heritage atmospheric chemistry models
-> tools
   -- It contains important tools to handle chemical file format, and visualization of chemical reaction networks
