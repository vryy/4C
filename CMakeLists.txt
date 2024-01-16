# Build system for baci and support tools.

#
### General Setup

cmake_minimum_required(VERSION 3.25)

# try to prevent modification of source directory
# note: some files may still be written before CMake can abort and need to be removed manually
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  message(
    FATAL_ERROR
      "In-source build not allowed. "
      "Please create a new directory, preferably next to the source directory, and run CMake from there. "
      "You may want to remove CMakeCache.txt and CMakeFiles/ which were created in the source directory."
    )
endif()

project(baci CXX C)

# Add flags that are vital to BACI's build system and cannot be disabled by users

# Export symbols (necessary for stacktraces)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -rdynamic")

# Print CMake version to screen
message(STATUS "Using CMake ${CMAKE_VERSION}")

# Enforce the C++ standard we are using and turn off compiler-specific extensions
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_CXX_EXTENSIONS Off)

# The version number.
set(baci_VERSION_MAJOR 1)
set(baci_VERSION_MINOR 0)

# specify where our own cmake modules are located
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)

# include all external (cmake) macros that we use
include(CheckCXXSourceRuns)
include(CheckCXXSourceCompiles)
include(FetchContent)

# include our own functions and macros
file(GLOB _function_files CONFIGURE_DEPENDS "cmake/functions/*.cmake")
foreach(_file ${_function_files})
  message(STATUS "Include ${_file}")
  include(${_file})
endforeach()

# If enabled, build all targets with address sanitizer
option(BACI_WITH_ADDRESS_SANITIZER "Compile BACI with address sanitizer" OFF)
if(BACI_WITH_ADDRESS_SANITIZER)
  add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
  add_link_options(-fsanitize=address -fno-omit-frame-pointer)
endif()

#
### Build in default values

if(NOT CMAKE_BUILD_TYPE)
  message(FATAL_ERROR "Choose a build type. Options are Debug and Release.")
endif()

###------------------------------------------------------------------- Default cmake settings
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

### process global options
include(cmake/setup_global_options.cmake)

# Run some system checks
file(GLOB _check_files CONFIGURE_DEPENDS "cmake/checks/*.cmake")
foreach(_file ${_check_files})
  message(STATUS "Include ${_file}")
  include(${_file})
endforeach()

###------------------------------------------------------------------ External Libraries

# Collect all external dependencies in a helper target
add_library(baci_all_enabled_external_dependencies INTERFACE)

baci_configure_dependency(HDF5 DEFAULT ON)
baci_configure_dependency(MPI DEFAULT ON)
baci_configure_dependency(Qhull DEFAULT ON)
baci_configure_dependency(Trilinos DEFAULT ON)
baci_configure_dependency(Boost DEFAULT ON)
baci_configure_dependency(ArborX DEFAULT OFF)
baci_configure_dependency(FFTW DEFAULT ON)
baci_configure_dependency(CLN DEFAULT ON)
baci_configure_dependency(MIRCO DEFAULT OFF)
baci_configure_dependency(Backtrace DEFAULT OFF)

# Generate the macro definition for all dependencies automatically
get_property(BACI_FLAGS_EXTERNAL_DEPENDENCIES GLOBAL PROPERTY BACI_FLAGS_EXTERNAL_DEPENDENCIES)
foreach(_dependency ${BACI_FLAGS_EXTERNAL_DEPENDENCIES})
  if(${_dependency})
    string(APPEND BACI_DEFINES_FOR_EXTERNAL_DEPENDENCIES "#define ${_dependency}\n")
  else()
    string(APPEND BACI_DEFINES_FOR_EXTERNAL_DEPENDENCIES "/* #undef ${_dependency} */\n")
  endif()
endforeach()

### end of configuration and feature detection

# Enable testing right away since we encounter src and tests interspersed in modules
enable_testing()

# NOTE: --mca orte_tmpdir_base [...] can be removed for openmpi >= 4.1.1 according to
# https://github.com/open-mpi/ompi/issues/8510#issuecomment-1329297350
set(MPIEXEC_EXTRA_OPTS_FOR_TESTING
    "--bind-to none --mca orte_tmpdir_base ${PROJECT_BINARY_DIR}/tmp"
    )

if(BACI_WITH_ADDRESS_SANITIZER)
  # Do not detect leaks, we only care about UB inducing memory issues
  string(APPEND MPIEXEC_EXTRA_OPTS_FOR_TESTING " -x LSAN_OPTIONS=detect_leaks=0")
  set(BACI_WITH_ADDRESS_SANITIZER_TEST_OPTIONS "LSAN_OPTIONS=detect_leaks=0")
endif()

set(FILTER_LIBRARIES post_common -Wl,--start-group ${BACI_LIBRARIES} -Wl,--end-group)

# hide libraries

if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/src)
endif(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
set(baciname baci-${BUILD_TYPE})

add_custom_target(post DEPENDS post_processor post_monitor)
add_custom_target(post_ensight DEPENDS post_processor)
add_custom_target(post_vtk DEPENDS post_processor)
add_custom_target(framework DEPENDS pre_exodus post)

set(FULL_TARGETS
    post_processor
    post_monitor
    pre_exodus
    pre_locsys
    create_rtd
    cut_test
    )

add_custom_target(full DEPENDS ${baciname} ${FULL_TARGETS})

include(cmake/setup_tests.cmake)

# Read in the source directory which defines all modules and executables
add_subdirectory(src)

# Process documentation
add_subdirectory(doc/doxygen EXCLUDE_FROM_ALL)
add_subdirectory(doc/readthedocs EXCLUDE_FROM_ALL)

###------------------------------------------------------------------ Tests

add_subdirectory(unittests)
add_subdirectory(tests/cut_test)
include(TestingFramework.cmake)

# Store the defined modules in a variable for consumption by baci_link_executable_with_all()
get_property(BACI_ALL_DEFINED_MODULE_TARGETS GLOBAL PROPERTY BACI_ALL_DEFINED_MODULE_TARGETS)

## After we defined all targets, we can now resolve the legacy linking
get_property(
  BACI_TARGETS_REQUIRING_CYCLIC_LEGACY_LINKING GLOBAL
  PROPERTY BACI_TARGETS_REQUIRING_CYCLIC_LEGACY_LINKING
  )
foreach(target ${BACI_TARGETS_REQUIRING_CYCLIC_LEGACY_LINKING})
  baci_link_executable_with_all(${target})
endforeach()
