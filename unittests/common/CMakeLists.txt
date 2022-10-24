# Define a special version of drt_error which will be linked to unit test executables
add_library(drt_error_for_testing dserror_for_testing.cpp)
target_link_libraries(drt_error_for_testing PUBLIC drt_lib_error_interface)
