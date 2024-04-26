# ==========================
# A collection of unit tests
# ==========================

configure_file(globals.cpp.in globals.cpp @ONLY)

# Enable testing functionality
enable_testing()

# Add test cases
setup_test(test_photolysis_box)

# copy yaml files
file(GLOB inputs *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/tests)
endforeach()
