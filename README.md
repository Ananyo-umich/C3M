Comprehensive Coupled Chemistry Model (C3M)

The code is developed by the Planetary Science Laboratory at the University of Michigan, Ann Arbor. It is an atmospheric photochemistry model that can be used as a standalone photochemistry model,
and can be coupled with hydrodynamic models for chemical transport in planetary environments.

Pre-requisite for running C3M:

-> Cantera
-> Athenapp
-> Eigen


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