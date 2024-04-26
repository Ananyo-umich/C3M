cmake_minimum_required(VERSION 3.18)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(
  C3M
  LANGUAGES CXX
  )

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message(STATUS "Default Build Type = Release")
  message(STATUS "")
else()
  message(STATUS "Build Type = ${CMAKE_BUILD_TYPE}")
  message(STATUS "")
endif()

option(PYTHON_BINDINGS "Build Python Bindings" ON)
if (PYTHON_BINDINGS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

# load all modules
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/)

## 1. set up compiler flags
message(STATUS "1. Set up project configuration ...")
message(STATUS "Include ${CMAKE_SOURCE_DIR}/cmake/compilers.cmake")
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/compilers.cmake)

## 2. set up project system libraries ##
message(STATUS "2. Set up project system libraries ...")

find_package(Eigen3 REQUIRED)
find_package(Cantera REQUIRED)
find_package(NetCDF REQUIRED)

## 3. set up project configure file and library ##
message(STATUS "3. Set up project libraries")
add_subdirectory(src)

## 4. set up examples and tests
message(STATUS "4. Set up unit tests")
#add_subdirectory(tests)

set(C3M_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL "C3M include directory")
