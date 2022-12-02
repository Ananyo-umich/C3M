
## set up all tests ##

file(GLOB _tests *.cpp)
foreach(_test ${_tests})
  foreach(build ${BUILD_TYPES})
    string(TOLOWER ${build} buildl)
    get_filename_component(_name ${_test} NAME_WE)

    add_executable(${_name}.${buildl} ${_test})

    set_target_properties(${_name}.${buildl}
      PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${EXECUTABLE_RELDIR}"
      COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${build}}
      )

    target_include_directories(${_name}.${buildl}
      SYSTEM PRIVATE
      ${C3M_INCLUDE_DIR}
      ${EIGEN3_INCLUDE_DIR}
      ${CANTERA_INCLUDE_DIR}
      ${ATHENAPP_INCLUDE_DIR}
      ${CMAKE_BINARY_DIR}
      )

    target_link_directories(${_name}.${buildl}
      SYSTEM PRIVATE
      ${CANTERA_LIBRARY_DIR}
      ${ATHENAPP_LIBRARY_DIR}
      )

    target_link_libraries(${_name}.${buildl}
      ${C3M_LIBRARY_${build}}
      ${CANTERA_LIBRARY}
      python3.9
      )
  endforeach()
endforeach()
