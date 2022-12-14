# Build system for baci and support tools.
# Including nightly test definitions.
#

#
### General Setup

cmake_minimum_required(VERSION 3.17)

project(baci CXX C Fortran)

# Print CMake version to screen
message("II Using CMake ${CMAKE_VERSION}")

# The version number.
set(baci_VERSION_MAJOR 1)
set(baci_VERSION_MINOR 0)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file(
  "${PROJECT_SOURCE_DIR}/src/headers/compile_settings.h.in"
  "${PROJECT_BINARY_DIR}/src/headers/compile_settings.h"
  )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}/src/headers")

# specify where our own cmake modules are located
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

# get current Baci revision from git
include(GetGitRevision)
get_git_revision_information()

# configure revision file to pass git revision information to the source code
configure_file(
  "${PROJECT_SOURCE_DIR}/src/headers/revision.H.in" "${PROJECT_BINARY_DIR}/src/headers/revision.H"
  )

# include all external (cmake) macros that we use
include(CheckCXXSourceRuns)
include(CheckCXXSourceCompiles)
include(FetchContent)

# include our own macros
include(BACILinkDefaultLibraries)
include(BACIAddLibrary)

#
### Build in default values

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE
      RelWithDebInfo
      CACHE STRING
            "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
            FORCE
      )
endif(NOT CMAKE_BUILD_TYPE)

#
### compile time options

if(${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
  add_definitions("-DDEBUG=ON")
endif()

option(
  DEBUGCUTLIBRARY "Turn on special output inside the cut library to be able to track errors" OFF
  )
if(DEBUGCUTLIBRARY)
  add_definitions("-DDEBUGCUTLIBRARY")
endif(DEBUGCUTLIBRARY)

option(
  TRILINOS_DEVELOP
  "Select Trilinos installation based on current develop branch (highly experimental!)"
  OFF
  )
if(TRILINOS_DEVELOP)
  message("II Trilinos installation: develop branch (highly experimental!)")
  add_definitions("-DTRILINOS_DEVELOP")
endif(TRILINOS_DEVELOP)

option(ENABLE_STACKTR "Enable printing of a stacktrace in case of a dserror" ON)
if(ENABLE_STACKTR)
  message("II Print stacktrace in dserror")
  add_definitions("-DENABLE_STACKTR")
endif(ENABLE_STACKTR)

# ignore selected cut test
option(IGNORE_UNSTABLE_CUTTEST "Ignore unstable cut test" OFF)
if(IGNORE_UNSTABLE_CUTTEST)
  message("-- Ignore unstable cut test ")
  add_definitions("-DIGNORE_UNSTABLE_CUTTEST")
endif(IGNORE_UNSTABLE_CUTTEST)

# fancy stuff

option(COLOROUTPUT "use colored output to the terminal" OFF)
if(COLOROUTPUT)
  add_definitions("-DCOLOROUTPUT")
endif(COLOROUTPUT)

option(DSERROR_DUMP "dserror creates a core file" OFF)
if(DSERROR_DUMP)
  add_definitions("-DDSERROR_DUMP")
endif(DSERROR_DUMP)

option(TRAP_FE "crash baci if a nan of inf occurs" ON)
if(TRAP_FE)
  add_definitions("-DTRAP_FE")
endif(TRAP_FE)

option(THROWELEMENTERRORS "throw errors at element level so all element error can be found" OFF)
if(THROWELEMENTERRORS)
  add_definitions("-DTHROWELEMENTERRORS")
endif(THROWELEMENTERRORS)

# Run some system checks
file(GLOB _check_files "cmake/checks/*.cmake")
foreach(_file ${_check_files})
  message(STATUS "Include ${_file}")
  include(${_file})
endforeach()

###------------------------------------------------------------------ External Libraries

find_package(
  HDF5
  COMPONENTS C HL
  REQUIRED
  )
find_package(MPI REQUIRED)
find_package(Qhull REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Trilinos REQUIRED)
find_package(
  Boost
  COMPONENTS graph system
  REQUIRED
  )
find_package(Netcdf REQUIRED)
find_package(ArborX)
find_package(FFTW)
find_package(ACML)
find_package(AMDLIBM)
find_package(CLN REQUIRED)

###------------------------------------------------------------------ List of files

### LEGACY
# Define variables that gather the names of all to-be-created libraries. These will be used
# for the executables that link against all libraries.

set(ELEMENT_LIBRARIES
    # cmake-format: sortable
    ele_ale
    ele_beam3
    ele_bele
    ele_discsh3
    ele_elemag
    ele_fluid
    ele_lubrication
    ele_membrane
    ele_porofluidmultiphase
    ele_rigidsphere
    ele_s8
    ele_scatra
    ele_so
    ele_so_plast
    ele_sp3
    ele_t3
    ele_t3cl
    ele_thermo
    ele_to3
    ele_w1
    )

set(MULTIFIELD_LIBRARIES
    # cmake-format: sortable
    alg_ehl
    alg_elch
    alg_fbi
    alg_fpsi
    alg_fs3i
    alg_fsi
    alg_fsi_xfem
    alg_immersed_problem
    alg_loma
    alg_pasi
    alg_poroelast
    alg_poromultiphase
    alg_poromultiphase_scatra
    alg_ssi
    alg_ssti
    alg_sti
    alg_tsi
    alg_two_phase
    alg_wear
    )

set(SINGLEFIELD_LIBRARIES
    # cmake-format: sortable
    alg_ale
    alg_artnet
    alg_contact
    alg_contact_aug
    alg_elemag
    alg_fluid
    alg_fluid_turbulence
    alg_lubrication
    alg_particle_algorithm
    alg_particle_engine
    alg_particle_interaction
    alg_particle_rigidbody
    alg_particle_wall
    alg_porofluidmultiphase
    alg_scatra
    alg_statinvanalysis
    alg_structure
    alg_structure_new
    alg_thermo
    alg_tutorial
    alg_xfluid
    alg_xstructure
    )

set(SUPPORT_LIBRARIES
    # cmake-format: sortable
    drt_fem_general
    drt_geometry_pair
    drt_inpar
    drt_io
    drt_levelset
    drt_lib
    drt_lib_artnet
    drt_lib_beamcontact
    drt_lib_beaminteraction
    drt_lib_binstrategy
    drt_lib_browniandyn
    drt_lib_constraint
    drt_lib_contact
    drt_lib_contact_constitutivelaw
    drt_lib_cut
    drt_lib_fiber
    drt_lib_geometric_search
    drt_lib_immersed
    drt_lib_io
    drt_lib_mortar
    drt_lib_poromulti_art_coupl
    drt_lib_statinv_ana
    drt_lib_struct_ale
    drt_lib_surfstress
    drt_lib_xfem
    drt_lib_xfem_manag
    drt_particleengine
    drt_patspec
    drt_red_airways
    )

set(BACI_LIBRARIES
    # cmake-format: sortable
    ${ELEMENT_LIBRARIES}
    ${MULTIFIELD_LIBRARIES}
    ${SINGLEFIELD_LIBRARIES}
    ${SUPPORT_LIBRARIES}
    drt_adapter
    drt_cardiovascular0d
    drt_comm
    drt_constraint
    drt_geometry
    drt_mat
    drt_matelast
    drt_mixture
    drt_mor
    drt_mortar
    drt_nurbs_discret
    drt_opti
    drt_stru_multi
    drt_timestepping
    drt_volmortar
    fortran
    headers
    linalg
    locacontinuation
    nlnsolver_nox
    pss_full
    rebalance
    solver
    )

set(FILTER_LIBRARIES post_drt_common -Wl,--start-group ${BACI_LIBRARIES} -Wl,--end-group)

###------------------------------------------------------------------ Binaries

# hide libraries

if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/src)
endif(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
set(baciname baci-${BUILD_TYPE})

add_custom_target(post DEPENDS post_processor post_drt_monitor)
add_custom_target(post_drt_ensight DEPENDS post_processor)
add_custom_target(post_drt_vtk DEPENDS post_processor)
add_custom_target(framework DEPENDS pre_exodus post)

set(FULL_TARGETS
    post_processor
    post_drt_monitor
    post_drt_s8convert
    pre_exodus
    pre_locsys
    cut_test
    )

add_custom_target(full DEPENDS ${baciname} ${FULL_TARGETS})

add_subdirectory(src)
add_subdirectory(doc EXCLUDE_FROM_ALL)

###------------------------------------------------------------------ Tests

enable_testing()
include(CTest)

###------------------------------------------------------------------ Unit Tests
add_subdirectory(tests/cut_test)
include(TestingFramework.cmake)
include(UnittestingFramework.cmake)
