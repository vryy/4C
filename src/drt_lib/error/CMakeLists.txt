# Define an interface library with the dserror declaration.
# This libray is necessary for the special implementation for testing.
add_library(drt_lib_error_interface INTERFACE)
target_include_directories(
  drt_lib_error_interface INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  )

# Define the production implementation.
add_library(drt_lib_error drt_dserror.cpp)
target_link_libraries(drt_lib_error PUBLIC drt_lib_error_interface)
