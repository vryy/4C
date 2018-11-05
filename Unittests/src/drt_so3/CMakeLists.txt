set_property(GLOBAL APPEND PROPERTY unit_files
${CMAKE_CURRENT_SOURCE_DIR}/unit_so_hex8.H
${CMAKE_CURRENT_SOURCE_DIR}/unit_so_tet4.H
${CMAKE_CURRENT_SOURCE_DIR}/unit_so_hex8_determinant_analysis.H
)

# Automatic version, but it does not detect new tests without rerunning cmake!
#file(GLOB LOCAL_UNIT_FILES "unit_*.H")
#set_property(GLOBAL APPEND PROPERTY unit_files ${LOCAL_UNIT_FILES})
