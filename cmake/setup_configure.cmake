## set up model configuration ##

if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
# set(CMAKE_BUILD_TYPE "DebugRelease")
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# float point operation precision
option(SinglePrecision "Enable single precision" OFF)

# equation of state
set(EquationOfState adiabatic_hydro
  CACHE STRING "Choose the equation of state for primitive-conserved conversion")
set_property(CACHE EquationOfState
  PROPERTY STRINGS
  adiabatic_hydro
  ideal_moist_hydro
  )

# coordinate system
set(CoordinateSystem cartesian
  CACHE STRING "Choose the coordinate system of the problem")
set_property(CACHE CoordinateSystem
  PROPERTY STRINGS
  cartesian
  cylindrical
  spherical_polar
  minkowski
  sinusoidal
  tilted
  schwarzschild
  kerr-schild
  gr_user
  )

# hydro flux solver
set(RiemannSolver hllc
  CACHE STRING "Choose the Riemann Solver")

# turbulence model
set(TurbulenceModel None
  CACHE STRING "Choose the turbulence model")
set_property(CACHE TurbulenceModel
  PROPERTY STRINGS
  None
  KEpsilon
  )

# ghost zone size
set(GhostZoneSize 3
  CACHE STRING "Set ghose zone size")

# array size
set(NumVapors 0
  CACHE STRING "Set number of vapors in the equation of state")
set(NumClouds 0
  CACHE STRING "Set number of clouds")
set(NumTracers 0
  CACHE STRING "Set number of fluid tracers")

# radiation stream
set(NumStreams 4
  CACHE STRING "Set number of radiation streams")
set(NumSpectralBands 0
  CACHE STRING "Set number of spectral bands")

# NetCDF output flag
option(UseNetCDF "Enable NetCDF output" ON)

# PNetCDF output flag
option(UsePNetCDF "Enable PNetCDF output" OFF)

# MPI flag
option(UseMPI "Enable MPI" OFF)

# CubedSphere flag
option(UseCubedSphere "Enable CubedSphere" OFF)

# HARP flag
option(UseHarp "Enable Harp package" OFF)

# CLIASTRO flag
option(UseCliAstro  "Enable CliAstro package" OFF)

# hydrostatic flag
option(Hydrostatic "Use hydrostatic log-pressure grid" OFF)
set(ReferencePressure 1.E5
  CACHE STRING "Set hydrostatic reference pressure")
set(PressureScaleHeight 30.E3
  CACHE STRING "Set hydrostatic pressure scale height")

set(ReferenceCloudPhase 2
  CACHE STRING "Set reference cloud phase for latent heat")

if (CMAKE_BUILD_TYPE MATCHES "Debug")
  if (NOT "DEBUG" IN_LIST BUILD_TYPES)
    list(APPEND BUILD_TYPES "DEBUG")
  endif()
endif()

if (CMAKE_BUILD_TYPE MATCHES "Release")
  if (NOT "RELEASE" IN_LIST BUILD_TYPES)
    list(APPEND BUILD_TYPES "RELEASE")
  endif()
endif()
