# Define a special version of drt_error which will be linked to unit test executables
add_library(drt_error_for_testing dserror_for_testing.cpp)
target_include_directories(drt_error_for_testing PRIVATE ${PROJECT_SOURCE_DIR})