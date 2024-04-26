# ==========================
# A collection of unit tests
# ==========================

configure_file(globals.cpp.in globals.cpp @ONLY)

# Enable testing functionality
enable_testing()

# Add test cases
setup_test(test_weno)
setup_test(test_glog)
setup_test(test_absorber)

setup_test(test_read_stellar_flux)
setup_test(test_read_cia_ff)
setup_test(test_read_rayleigh)
setup_test(test_read_yaml)
setup_test(test_helios_ck_read)
setup_test(test_cs_velocity_rotation)

# copy yaml files
file(GLOB inputs *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/tests)
endforeach()
