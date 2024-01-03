# Kokkos is typically pulled in via Trilinos. If no location has been given,
# try the same location as Trilinos. If no Trilinos location exists, users
# will get an error to provide that one first.
if(Trilinos_ROOT AND NOT Kokkos_ROOT)
  set(Kokkos_ROOT
      ${Trilinos_ROOT}
      CACHE PATH "Path to Kokkos installation"
      )
endif()

# We only support Trilinos versions that provide a config file.
find_package(Trilinos REQUIRED)

message(STATUS "Trilinos packages: ${Trilinos_PACKAGE_LIST}")
message(STATUS "Trilinos TPLs: ${Trilinos_TPL_LIST}")

baci_add_dependency(baci_all_enabled_external_dependencies Trilinos::all_selected_libs)
