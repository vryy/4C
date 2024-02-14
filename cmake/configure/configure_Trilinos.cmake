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

if(Trilinos_FOUND AND NOT TARGET Trilinos::all_selected_libs)
  # In preparation for newer Trilinos releases, create a target
  # Trilinos::all_selected_libs with the correct dependencies
  add_library(Trilinos::all_selected_libs IMPORTED INTERFACE)
  target_include_directories(
    Trilinos::all_selected_libs
    SYSTEM
    INTERFACE ${Trilinos_INCLUDE_DIRS}
    INTERFACE ${Trilinos_TPL_INCLUDE_DIRS}
    )
  target_link_libraries(
    Trilinos::all_selected_libs INTERFACE ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
    )
endif()

target_link_libraries(baci_all_enabled_external_dependencies INTERFACE Trilinos::all_selected_libs)
