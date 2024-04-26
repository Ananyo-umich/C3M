# compile c3m

set(namel c3m)
string(TOUPPER ${namel} nameu)

if (PYTHON_BINDINGS)
  set(compile_flags_ "${${PROJECT_NAME}_C_FLAGS} -fPIC")
else()
  set(compile_flags_ "${${PROJECT_NAME}_C_FLAGS}")
endif()

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
  COMPILE_FLAGS ${compile_flags_}
  )

target_include_directories(${namel}_${buildl}
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${EIGEN3_INCLUDE_DIR}
  ${CANTERA_INCLUDE_DIR}
  )

target_link_libraries(${namel}_${buildl}
  ${CANTERA_LIBRARY}
  )

add_library(c3m::c3m ALIAS ${namel}_${buildl})

set(C3M_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE INTERNAL "C3M include directory")
