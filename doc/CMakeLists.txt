#
# Build routine/test for doxygen documentation

#
### General Setup

cmake_minimum_required(VERSION 2.8)
project(documentation CXX C)

# generates documentation in build folder
# Dependencies: Doxyfile <- Doxyfile.defs <- gen_doxygen_defs.py
add_custom_command ( OUTPUT Doxyfile COMMAND cat ${CMAKE_SOURCE_DIR}/doc/Doxyfile.tmpl Doxyfile.defs > Doxyfile && sed -i -e "s/####REVISION####/${baci_WC_REVISION}/" -e "s:####SOURCE####:${CMAKE_SOURCE_DIR}:" Doxyfile MAIN_DEPENDENCY ${PROJECT_BINARY_DIR}/Doxyfile.defs )
add_custom_command ( OUTPUT Doxyfile.defs COMMAND ${CMAKE_SOURCE_DIR}/gen_doxygen_defs.py "${PROJECT_BINARY_DIR}/../" "${PROJECT_BINARY_DIR}" )


# start build with "make documentation"
# N.B. this target is always considered out of date and therefore doxygen always runs
add_custom_target( documentation ALL echo COMMAND doxygen Doxyfile DEPENDS ${PROJECT_BINARY_DIR}/Doxyfile )


