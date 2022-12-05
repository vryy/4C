add_library(fortran INTERFACE)
target_include_directories(fortran INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
