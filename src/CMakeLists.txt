# compile c3m

set(namel c3m)
string(TOUPPER ${namel} nameu)

file(GLOB _src_files
  *.cpp
  )

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

add_library(${namel}_${buildl}
  STATIC
  ${_src_files}
  )

set_target_properties(${namel}_${buildl}
  PROPERTIES
  COMPILE_FLAGS ${C3M_CXX_FLAGS}
  )

target_include_directories(${namel}_${buildl}
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${EIGEN3_INCLUDE_DIR}
  ${CANTERA_INCLUDE_DIR}
  ${HDF5_INCLUDE_DIR}
  ${NETCDF_INCLUDES}
  )

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  set(extra_flag stdc++fs)
endif()

target_link_libraries(${namel}_${buildl}
  ${CANTERA_LIBRARY}
  ${NETCDF_LIBRARIES}
  ${CMAKE_DL_LIBS}
  ${HDF5_LIBRARIES}
  ${extra_flag}
  )

add_library(c3m::c3m ALIAS ${namel}_${buildl})
