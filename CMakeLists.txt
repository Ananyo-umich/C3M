## 0.set up project ##

cmake_minimum_required(VERSION 3.18)
message(STATUS "Set up canoe project ...")

project(
  C3M
  VERSION 0.0.1
  LANGUAGES C CXX Fortran
  )
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)

# 1. set up compiler flags
message(STATUS "")
message(STATUS "Set up project configuration ...")
message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags.cmake")
include(${CMAKE_SOURCE_DIR}/cmake/setup_compiler_flags.cmake)

## 2. set up project configuration ##
message(STATUS "")
message(STATUS "Set up project configuration ...")
message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/setup_configure.cmake")
include(${CMAKE_SOURCE_DIR}/cmake/setup_configure.cmake)

## 3. set up project paths ##
message(STATUS "")
message(STATUS "Set up project paths ...")
set(EXECUTABLE_RELDIR "bin")
set(INCLUDE_RELDIR "include")
set(SOURCE_RELDIR "source")

## 4. set up project system libraries ##
message(STATUS "")
message(STATUS "Set up project system libraries ...")

find_package(Eigen3 3.3 REQUIRED)
find_package(Cantera REQUIRED)
find_package(HDF5 REQUIRED)
find_package(OpenMP REQUIRED)  # Find the OpenMP package
if (${UseNetCDF})
  find_package(NetCDF REQUIRED)
endif()
if (${UsePNetCDF})
  find_package(PNetCDF REQUIRED)
endif()

if (${UseMPI})
  find_package(MPI REQUIRED)
  set(MPIOption MPI_PARALLEL)
else()
  set(MPIOption NOT_MPI_PARALLEL)
endif()

if (${UseCubedSphere})
  set(CubedSphereOption CUBED_SPHERE)
else()
  set(CubedSphereOption NOT_CUBED_SPHERE)
endif()

## 5. set up project external libraries ##
message(STATUS "")
message(STATUS "Set up project external libraries ...")

# compile athenapp
if (${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
  add_subdirectory(athenapp)
endif()

## 6. set up compilation and installation ##
message(STATUS "")
message(STATUS "Set up project complication and installation ...")

# configure file
configure_file(${CMAKE_SOURCE_DIR}/configure.hpp.in configure.hpp @ONLY)

#add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(tests)
