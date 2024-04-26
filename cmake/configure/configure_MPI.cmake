# This find module is provided by CMake

# Disable deprecated CXX bindings in MPI. These often lead to compiler warnings.
set(MPI_CXX_SKIP_MPICXX ON)
find_package(MPI REQUIRED)

target_link_libraries(four_c_all_enabled_external_dependencies INTERFACE MPI::MPI_CXX)
