# ==========================
# A collection of unit tests
# ==========================

configure_file(globals.cpp.in globals.cpp @ONLY)

# Enable testing functionality
enable_testing()

# Add test cases
setup_test(test_photolysis_box)
setup_test(test_photochem_1d)
setup_test(test_advection_diffusion)
setup_test(test_chapman)
setup_test(test_SOx_network)
setup_test(test_Jupiter_diffusion)
setup_test(PhotoChemBox)
setup_test(PhotoChem1D)

# copy yaml files
file(GLOB inputs *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/tests)
endforeach()
