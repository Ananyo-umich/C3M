# =====================
# Examples for C3M wip
# =====================


# 1. Compile the problem
setup_problem(PhotoChemBox)

# 2. Copy input files to run directory
file(GLOB inputs *.yaml)
foreach(input ${inputs})
    file(COPY ${input} DESTINATION ${CMAKE_BINARY_DIR}/bin)
endforeach()
