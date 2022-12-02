# Installs C3M library
#
#   library: c3m_debug
#   library: c3m_release
#
# Define the following variables
#
#   C3M_INCLUDE_DIR
#   C3M_LIBRARY_DEBUG
#   C3M_LIBRARY_RELEASE
#

file(GLOB _src_files
  *.cpp
  )

foreach(build ${BUILD_TYPES})
  string(TOLOWER ${build} buildl)

  add_library(c3m_${buildl} ${_src_files})

  set_target_properties(c3m_${buildl}
    PROPERTIES
    COMPILE_FLAGS ${CMAKE_CXX_FLAGS_${build}}
    )

  target_include_directories(c3m_${buildl}
    SYSTEM PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${CANTERA_INCLUDE_DIR}
    ${CMAKE_BINARY_DIR}
    )

  target_link_directories(c3m_${buildl}
    SYSTEM PRIVATE
    ${CANTERA_LIBRARY_DIR}
    )

  target_link_libraries(c3m_${buildl}
    ${CANTERA_LIBRARY}
    python3.9
    )

  set(C3M_LIBRARY_${build} c3m_${buildl}
    CACHE STRING 
    "c3m library ${build}")
  mark_as_advanced(C3M_LIBRARY_${build})
endforeach()

set(C3M_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  CACHE PATH 
  "C3M include directory")
mark_as_advanced(C3M_INCLUDE_DIR)
