# A small macro used for setting up the build of a problem.
#
# Usage: setup_problem(name)

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

macro(setup_problem namel)
  add_executable(${namel}.${buildl} ${namel}.cpp)
  
  set_target_properties(${namel}.${buildl} PROPERTIES COMPILE_FLAGS
                                                      ${C3M_CXX_FLAGS})

  
  target_include_directories(
    ${namel}.${buildl} PRIVATE ${CMAKE_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR} ${MPI_CXX_INCLUDE_PATH} ${NETCDF_INCLUDES} ${PNETCDF_INCLUDE_DIR} ${CANTERA_INCLUDE_DIR} ${YAMLPP_INCLUDE_DIR})

  target_link_libraries(${namel}.${buildl}
	  PRIVATE ${C3M_LIBRARIES} ${EIGEN3_LIBRARIES}
	  ${MPI_CXX_LIBRARIES} ${NETCDF_LIBRARIES} ${PNETCDF_LIBRARIES} ${CANTERA_LIBRARIES} ${YAMLPP_LIBRARIES} )


endmacro()
