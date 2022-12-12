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
include(Deprecated_BACIAddLibraryLinkedWithEverything)
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
set(OBJS_POST_DRT_MONITOR src/post_monitor/post_drt_monitor.cpp)
set(OBJS_DRT_UTILS
    src/drt_fluid/fluid_functions.cpp
    src/drt_fluid/fluid_rotsym_periodicbc_utils.cpp
    src/drt_fluid_xfluid/xfluid_functions.cpp
    src/drt_fluid_xfluid/xfluid_functions_combust.cpp
    src/drt_lib/drt_condition_selector.cpp
    src/drt_lib/drt_condition_utils.cpp
    src/drt_lib/drt_dirichletextractor.cpp
    src/drt_lib/drt_elementdefinition.cpp
    src/drt_lib/drt_elementtype.cpp
    src/drt_lib/drt_exporter.cpp
    src/drt_lib/drt_matchingoctree.cpp
    src/drt_lib/drt_parser.cpp
    src/drt_lib/drt_resulttest.cpp
    src/drt_lib/drt_utils.cpp
    src/drt_lib/drt_utils_parallel.cpp
    src/drt_lib/epetra_utils.cpp
    src/drt_beaminteraction/drt_utils_parallel_proctoproc.cpp
    src/rebalance/drt_utils_rebalancing.cpp
    src/drt_structure_new/str_functions.cpp
    src/drt_lib/voigt_notation.cpp
    )
set(OBJS_PRE_EXODUS
    src/pre_exodus/pre_exodus.cpp
    src/pre_exodus/pre_exodus_centerline.cpp
    src/pre_exodus/pre_exodus_readbc.cpp
    src/pre_exodus/pre_exodus_reader.cpp
    src/pre_exodus/pre_exodus_soshextrusion.cpp
    src/pre_exodus/pre_exodus_validate.cpp
    src/pre_exodus/pre_exodus_writedat.cpp
    )
set(OBJS_POST_DRT_COMMON src/post_drt_common/post_drt_common.cpp)
set(OBJS_POST_PROCESSOR_MAIN
    src/post_drt_common/post_processor.cpp
    src/post_drt_common/post_single_field_writers.cpp
    src/post_drt_common/post_structure_stress.cpp
    src/post_drt_common/post_thermo_heatflux.cpp
    src/post_drt_common/post_filter_base.cpp
    src/post_drt_common/post_writer_base.cpp
    )
set(OBJS_POST_DRT_ENSIGHT
    src/post_ensight/post_drt_ensight_writer.cpp src/post_ensight/post_drt_ensight_writer_nurbs.cpp
    )
set(OBJS_POST_DRT_VTK
    src/post_vtk/post_drt_vtu_writer.cpp
    src/post_vtk/post_drt_vtu_writer_node_based.cpp
    src/post_vtk/post_drt_vti_writer.cpp
    src/post_vtk/post_drt_vtk_writer.cpp
    )
set(OBJS_BACI_MAIN
    src/global_full/global_init_control.cpp
    src/global_full/global_inp_control.cpp
    src/global_full/global_cal_control.cpp
    src/global_full/global_control.cpp
    src/global_full/baci.cpp
    )
set(OBJS_PRE_LOCSYS src/pre_locsys/pre_locsys.cpp)
set(OBJS_POST_DRT_S8CONVERT src/post_s8convert/post_drt_s8convert.cpp)
set(OBJS_CUT_TEST
    tests/cut_test/cut_test_double.cpp
    tests/cut_test/cut_test_geometry.cpp
    tests/cut_test/cut_test_intersection.cpp
    tests/cut_test/cut_test_levelset.cpp
    tests/cut_test/cut_test_levelset_no_loc_coord.cpp
    tests/cut_test/cut_test_loader.cpp
    tests/cut_test/cut_test_main.cpp
    tests/cut_test/cut_test_selfcut.cpp
    tests/cut_test/cut_test_simple.cpp
    tests/cut_test/cut_test_triangulateFacet.cpp
    tests/cut_test/cut_test_units.cpp
    tests/cut_test/cut_test_utils.cpp
    tests/cut_test/cut_test_fluidfluid.cpp
    tests/cut_test/cut_test_newlineimplementation.cpp
    tests/cut_test/cut_test_alex_39.cpp
    tests/cut_test/cut_test_alex_40.cpp
    tests/cut_test/cut_test_alex_41.cpp
    tests/cut_test/cut_test_alex_42.cpp
    tests/cut_test/cut_test_alex_43.cpp
    tests/cut_test/cut_test_alex_44.cpp
    tests/cut_test/cut_test_alex_45.cpp
    tests/cut_test/cut_test_alex_46.cpp
    tests/cut_test/cut_test_alex_47.cpp
    tests/cut_test/cut_test_alex_48.cpp
    tests/cut_test/cut_test_alex_49.cpp
    tests/cut_test/cut_test_alex_50.cpp
    tests/cut_test/cut_test_alex_51.cpp
    tests/cut_test/cut_test_alex_52.cpp
    tests/cut_test/cut_test_alex_53.cpp
    tests/cut_test/cut_test_alex_54.cpp
    tests/cut_test/cut_test_alex_55.cpp
    tests/cut_test/cut_test_alex_56.cpp
    tests/cut_test/cut_test_alex_57.cpp
    tests/cut_test/cut_test_alex_58.cpp
    tests/cut_test/cut_test_alex_59.cpp
    tests/cut_test/cut_test_alex_60.cpp
    tests/cut_test/cut_test_alex_61.cpp
    tests/cut_test/cut_test_alex_62.cpp
    tests/cut_test/cut_test_axel_8.cpp
    tests/cut_test/cut_test_axel_9.cpp
    tests/cut_test/cut_test_axel_10.cpp
    tests/cut_test/cut_test_benedikt_1.cpp
    tests/cut_test/cut_test_hex8_twintri.cpp
    tests/cut_test/cut_test_shadan_6.cpp
    tests/cut_test/cut_test_volume.cpp
    tests/cut_test/drt_dserror.cpp
    tests/cut_test/cut_test_intersection.cpp
    tests/cut_test/cut_test_bacigenerated_1890.cpp
    tests/cut_test/cut_test_bacigenerated_6890.cpp
    tests/cut_test/cut_test_bacigenerated_1858.cpp
    tests/cut_test/cut_test_bacigenerated_1970.cpp
    tests/cut_test/cut_test_bacigenerated_1910.cpp
    tests/cut_test/cut_test_bacigenerated_2010.cpp
    tests/cut_test/cut_test_bacigenerated_1860.cpp
    tests/cut_test/cut_test_facet_failed_cln.cpp
    tests/cut_test/cut_test_bacigenerated_6920.cpp
    tests/cut_test/cut_test_bacigenerated_7019.cpp
    tests/cut_test/cut_test_bacigenerated_6923.cpp
    tests/cut_test/cut_test_bacigenerated_41534.cpp
    tests/cut_test/cut_test_bacigenerated_43244.cpp
    tests/cut_test/cut_test_bacigenerated_463638.cpp
    tests/cut_test/cut_test_bacigenerated_622320.cpp
    tests/cut_test/cut_test_bacigenerated_227469.cpp
    tests/cut_test/cut_test_bacigenerated_622829.cpp
    tests/cut_test/cut_test_bacigenerated_627558.cpp
    tests/cut_test/cut_test_bacigenerated_7022.cpp
    tests/cut_test/cut_test_bacigenerated_238425.cpp
    tests/cut_test/cut_test_bacigenerated_197489.cpp
    tests/cut_test/cut_test_bacigenerated_79216.cpp
    tests/cut_test/cut_test_bacigenerated_369096.cpp
    tests/cut_test/cut_test_bacigenerated_238343.cpp
    tests/cut_test/cut_test_bacigenerated_26182.cpp
    )
set(OBJS_GRAPH_TEST
    tests/graph_test/graph_test_main.cpp
    tests/graph_test/graph_test_read.cpp
    tests/graph_test/graph_test_gnuplot.cpp
    tests/graph_test/graph_test_output.cpp
    tests/graph_test/graph_test_find_cycles.cpp
    )
### LEGACY
# Define variables that gather the names of all to-be-created libraries. These will be used within
# deprecated_baci_add_library_linked_with_everything() to link every target with everything.

set(ELEMENT_LIBRARIES
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
    ele_sp3
    ele_t3
    ele_t3cl
    ele_thermo
    ele_to3
    ele_w1
    )

set(MULTIFIELD_LIBRARIES
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
    alg_ale
    alg_artnet
    alg_contact
    alg_contact_aug
    alg_elemag
    alg_fluid
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
    # Libraries needed by elements
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
    ${ELEMENT_LIBRARIES}
    ${MULTIFIELD_LIBRARIES}
    ${SINGLEFIELD_LIBRARIES}
    ${SUPPORT_LIBRARIES}
    drt_adapter
    drt_cardiovascular0d
    drt_comm
    drt_constraint
    drt_geometry
    drt_lib_register
    drt_lib_utils
    drt_mat
    drt_matelast
    drt_mixture
    drt_mor
    drt_mortar
    drt_nurbs_discret
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

###------------------------------------------------------------------ Binaries

# hide libraries

if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/src)
endif(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

add_library(headers INTERFACE)
target_include_directories(headers INTERFACE ${PROJECT_SOURCE_DIR}/src/headers)

deprecated_baci_add_library_linked_with_everything(drt_lib_utils ${OBJS_DRT_UTILS})

# filter framework
deprecated_baci_add_library_linked_with_everything(post_drt_common ${OBJS_POST_DRT_COMMON})
deprecated_baci_add_library_linked_with_everything(post_drt_ensight_lib ${OBJS_POST_DRT_ENSIGHT})
deprecated_baci_add_library_linked_with_everything(post_drt_vtk_lib ${OBJS_POST_DRT_VTK})

target_link_libraries(post_drt_vtk_lib PUBLIC post_drt_common)
target_link_libraries(post_drt_ensight_lib PUBLIC post_drt_common)

string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
set(baciname baci-${BUILD_TYPE})

# the baci executable
add_executable(${baciname} ${OBJS_BACI_MAIN})

add_executable(post_processor EXCLUDE_FROM_ALL ${OBJS_POST_PROCESSOR_MAIN})
add_executable(post_drt_monitor EXCLUDE_FROM_ALL ${OBJS_POST_DRT_MONITOR})
add_executable(post_drt_s8convert EXCLUDE_FROM_ALL ${OBJS_POST_DRT_S8CONVERT})
add_executable(pre_exodus EXCLUDE_FROM_ALL ${OBJS_PRE_EXODUS})
add_executable(pre_locsys EXCLUDE_FROM_ALL ${OBJS_PRE_LOCSYS})

target_link_libraries(${baciname} ${BACI_LIBRARIES})

set(FILTER_LIBRARIES post_drt_common ${BACI_LIBRARIES})

target_link_libraries(post_processor post_drt_ensight_lib post_drt_vtk_lib ${FILTER_LIBRARIES})
target_link_libraries(post_drt_monitor ${FILTER_LIBRARIES})
target_link_libraries(post_drt_s8convert ${FILTER_LIBRARIES})
target_link_libraries(pre_locsys ${FILTER_LIBRARIES})
target_link_libraries(pre_exodus ${FILTER_LIBRARIES})

add_custom_target(baci DEPENDS ${baciname})
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

# various post-processing tools
add_custom_command(
  TARGET post_processor
  POST_BUILD
  COMMAND
    ${CMAKE_COMMAND} -D CMAKE_PROJECT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR} -P
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/scripts/baci_create_post_scripts.cmake
  )

add_subdirectory(src)

add_subdirectory(doc EXCLUDE_FROM_ALL)

###------------------------------------------------------------------ Tests

enable_testing()
include(CTest)

###------------------------------------------------------------------ Unit Tests
add_executable(cut_test EXCLUDE_FROM_ALL ${OBJS_CUT_TEST})
target_link_libraries(
  cut_test
  drt_lib_cut
  ${QHULL_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BACI_LIBRARIES}
  )

include(TestingFramework.cmake)
include(UnittestingFramework.cmake)
