# A small macro used for setting up the build of a test.
#
# Usage: setup_test(name)

macro(setup_test namel)
  string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
  string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

  add_executable(${namel}.${buildl} ${namel}.cpp globals.cpp)

  set_target_properties(${namel}.${buildl} PROPERTIES COMPILE_FLAGS
                                                      ${C3M_CXX_FLAGS})

  target_include_directories(${namel}.${buildl} PRIVATE ${CMAKE_SOURCE_DIR}
                                                        ${CANTERA_INCLUDE_DIR})

  target_link_libraries(
    ${namel}.${buildl} PRIVATE gtest_main ${CANTERA_LIBRARY}
                               application_${buildl} c3m::c3m)

  add_test(NAME ${namel}.${buildl} COMMAND ${namel}.${buildl})
endmacro()
