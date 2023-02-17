# Define an interface library with the dserror declaration.
# This libray is necessary for the special implementation for testing.
add_library(lib_error_interface INTERFACE)
target_include_directories(
  lib_error_interface INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

# Define the production implementation.
add_library(lib_error lib_dserror.cpp)

set(LIBS_LIB_ERROR
    # cmake-format: sortable
    lib_error_interface
    )

baci_add_dependency(lib_error ${LIBS_LIB_ERROR})
