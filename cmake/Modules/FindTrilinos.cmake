# temporarily disable our own finders for module lookup
list(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

# Get Trilinos as one entity
set(CMAKE_PREFIX_PATH ${Trilinos_PREFIX})
set(Trilinos_DIR "${Trilinos_PREFIX}/lib/cmake/Trilinos")

# call the built-in finder with the specifications set outside
find_package(Trilinos)

# re-enable our own finders
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

# Echo trilinos build info just for fun
message(STATUS "Found Trilinos: ${Trilinos_DIR} (Version ${Trilinos_VERSION})")
message(STATUS "Trilinos packages: ${Trilinos_PACKAGE_LIST}")
message(STATUS "Trilinos TPLs: ${Trilinos_TPL_LIST}")

# get trilinos version information
include(GetTrilinosVersion)
get_trilinos_version()

# configure trilinos version file to pass information to the source code
configure_file(
  "${PROJECT_SOURCE_DIR}/src/headers/trilinos_version.H.in"
  "${PROJECT_BINARY_DIR}/src/headers/trilinos_version.H"
  )

if(Trilinos_FOUND AND NOT TARGET Trilinos::all_selected_libs)
  # In preparation for newer Trilinos releases, create a target
  # Trilinos::all_selected_libs with the correct dependencies
  add_library(_Trilinos_all_selected_libs INTERFACE)
  target_include_directories(
    _Trilinos_all_selected_libs
    SYSTEM
    INTERFACE ${Trilinos_INCLUDE_DIRS}
    INTERFACE ${Trilinos_TPL_INCLUDE_DIRS}
    )
  target_link_libraries(
    _Trilinos_all_selected_libs INTERFACE ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES}
    )

  add_library(Trilinos::all_selected_libs ALIAS _Trilinos_all_selected_libs)
endif()

list(APPEND BACI_ALL_ENABLED_EXTERNAL_LIBS Trilinos::all_selected_libs)
