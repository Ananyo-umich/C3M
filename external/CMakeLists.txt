###############################
## set up external libraries ##
###############################

## extract CARMA
#
if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/CARMA30/CMakeLists.txt)
  message(STATUS "Found CARMA ...")
else()
  message(STATUS "Extracting CARMA...")
  file(ARCHIVE_EXTRACT
    INPUT ${CMAKE_CURRENT_SOURCE_DIR}/CARMA3.0.tar.gz
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})
endif()
add_subdirectory(CARMA30)
