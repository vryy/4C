# get trilinos prefixes for available packages
function(add_trilinos_package_definitions target trilinos_package_list)
  foreach(package ${trilinos_package_list})
    string(TOUPPER ${package} PACKAGE)
    set(PREFIX "BACI_TRILINOS_WITH_${PACKAGE}")
    target_compile_definitions(${target} INTERFACE ${PREFIX}=1)
  endforeach()
endfunction()

# temporarily disable our own finders for module lookup
list(REMOVE_ITEM CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# Get Trilinos as one entity
set(CMAKE_PREFIX_PATH ${Trilinos_PREFIX})
set(Trilinos_DIR "${Trilinos_PREFIX}/lib/cmake/Trilinos")

# call the built-in finder with the specifications set outside
find_package(Trilinos)

# re-enable our own finders
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# Echo trilinos build info just for fun
message(STATUS "Found Trilinos: ${Trilinos_DIR} (Version ${Trilinos_VERSION})")
message(STATUS "Trilinos packages: ${Trilinos_PACKAGE_LIST}")
message(STATUS "Trilinos TPLs: ${Trilinos_TPL_LIST}")

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
  add_trilinos_package_definitions(_Trilinos_all_selected_libs "${Trilinos_PACKAGE_LIST}")
  add_library(Trilinos::all_selected_libs ALIAS _Trilinos_all_selected_libs)
else()
  add_trilinos_package_definitions(Trilinos::all_selected_libs "${Trilinos_PACKAGE_LIST}")
endif()

list(APPEND BACI_ALL_ENABLED_EXTERNAL_LIBS Trilinos::all_selected_libs)
