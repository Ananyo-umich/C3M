
# set up project

cmake_minimum_required(VERSION 3.18)

message(STATUS "Set up compiler ...")

project(
  C3M
  VERSION 0.0.1
  LANGUAGES C CXX Fortran
  )
set(CMAKE_CXX_STANDARD 14)

message(STATUS "Done.")

message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/setup_configure.cmake")
include(${CMAKE_SOURCE_DIR}/cmake/setup_configure.cmake)

# set up model parameters
message(STATUS "")
message(STATUS "Set up project configuration ...")

# set up configure file
configure_file(configure.hpp.in configure.hpp @ONLY)

message(STATUS "Done.")

## set up project paths ##

message(STATUS "")
message(STATUS "Set up project paths ...")

set(EXECUTABLE_RELDIR "bin")
set(INCLUDE_RELDIR "include")
set(LIBRARY_RELDIR "lib${PACKAGE_NAME}")
set(SHARE_RELDIR "share/${PACKAGE_NAME}")

message(STATUS "Done.")


## debug and/or release mode ##
message(STATUS "")
message(STATUS "Set up compiler flags ...")
include(cmake/setup_compiler_flags.cmake)

if (CMAKE_BUILD_TYPE MATCHES "Debug")
  list(APPEND BUILD_TYPES "DEBUG")
endif()

if (CMAKE_BUILD_TYPE MATCHES "Release")
  list(APPEND BUILD_TYPES "RELEASE")
endif()

foreach(build ${BUILD_TYPES})
  message(STATUS "build = ${build}")
  message(STATUS "C compiler flags = ${CMAKE_C_FLAGS_${build}}")
  message(STATUS "C++ compiler flags = ${CMAKE_CXX_FLAGS_${build}}")
  message(STATUS "Fortran compiler flags = ${CMAKE_Fortran_FLAGS_${build}}")
endforeach()

message(STATUS "Done.")

## set up system libraries ##
# use cantera
set(CANTERA_INCLUDE_DIR /home/ananyo/opt/include)
set(CANTERA_LIBRARY_DIR /home/ananyo/opt/lib)
set(CANTERA_LIBRARY cantera)

## compilation and installation ##

message(STATUS "")
message(STATUS "Set up complication and installation ...")

add_subdirectory(external)
add_subdirectory(src)
add_subdirectory(tests)

message(STATUS "Done.")
message(STATUS "")
