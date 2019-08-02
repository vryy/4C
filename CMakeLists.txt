#
# Build system for baci and support tools.
# Including nightly test definitions.
#

#
### General Setup

cmake_minimum_required(VERSION 2.8)
project(baci CXX C Fortran)

# The version number.
set (baci_VERSION_MAJOR 1)
set (baci_VERSION_MINOR 0)

# configure a header file to pass some of the CMake settings
# # to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/src/headers/compile_settings.h.in"
  "${PROJECT_BINARY_DIR}/src/headers/compile_settings.h"
  )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}/src/headers")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# get current Baci revision from git
include(GetGitRevision)
get_git_revision_information()

# configure revision file to pass git revision information to the source code
configure_file(
  "${PROJECT_SOURCE_DIR}/src/headers/revision.H.in"
  "${PROJECT_BINARY_DIR}/src/headers/revision.H")

# include all external (cmake) macros that we use
INCLUDE(CheckCXXSourceRuns)
INCLUDE(CheckCXXSourceCompiles)

# define macro to remove duplicate libraries from the link line
MACRO(REMOVE_DUPLICATES _list)
  IF(NOT "${${_list}}" STREQUAL "")
    LIST(REVERSE ${_list})
    LIST(REMOVE_DUPLICATES ${_list})
    LIST(REVERSE ${_list})
  ENDIF()
ENDMACRO()

#
### Build in default values

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#
# Legacy defines

add_definitions("-DLINUX_MUENCH")
add_definitions("-DPARALLEL") # TODO remove PARALLEL flag? exact meaning of this flag?
add_definitions("-DHAVE_PARMETIS") # TODO remove PARMETIS flag
add_definitions("-DD_SHELL8") # TODO remove D_SHELL8 flag

# external library variable
set(LIBRARIES "")

#
### compile time options

if (${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
  add_definitions("-DDEBUG=ON")
endif ()

option (DEBUGCUTLIBRARY "Turn on special output inside the cut library to be able to track errors" OFF)
if (DEBUGCUTLIBRARY)
  add_definitions("-DDEBUGCUTLIBRARY")
endif (DEBUGCUTLIBRARY)

option (HAVE_TEKO "Turn on Teko support (needs Teko and Stratimikos to be compiled within Trilinos)" OFF)
if (HAVE_TEKO)
  MESSAGE("II Turn on support for Teko")
  add_definitions("-DHAVE_TEKO")
endif (HAVE_TEKO)

option (HAVE_MueLu "Turn on MueLu support (needs MueLu to be compiled within Trilinos" OFF)
if (HAVE_MueLu)
  MESSAGE("II Turn on support for MueLu")
  add_definitions("-DHAVE_MueLu")
endif (HAVE_MueLu)

option (HAVE_Stratimikos "Turn on Stratimikos support (needs Stratimikos to be compiled within Trilinos" OFF)
if (HAVE_Stratimikos)
  MESSAGE("II Turn on support for Stratimikos")
  add_definitions("-DHAVE_Stratimikos")
endif (HAVE_Stratimikos)

option (HAVE_Kokkos "Turn on Kokkos support (needs Kokkos to be compiled with Trilinos)" OFF)
if (HAVE_Kokkos)
  MESSAGE("II Turn on support for Kokkos")
  add_definitions("-DHAVE_Kokkos")
endif(HAVE_Kokkos)

option (ENABLE_STACKTR "Enable printing of a stacktrace in case of a dserror" ON)
if (ENABLE_STACKTR)
  MESSAGE("II Print stacktrace in dserror")
  add_definitions("-DENABLE_STACKTR")
endif (ENABLE_STACKTR)


#   additional features


find_package(HDF5 REQUIRED)
if (HDF5_FOUND)
  MESSAGE("-- HDF5 include directory: ${HDF5_INCLUDE_DIRS}")
  MESSAGE("-- HDF5 library directory: ${HDF5_LIBRARY}")
  MESSAGE("-- HDF5 libraries: ${HDF5_LIBRARIES}")
  include_directories(${HDF5_INCLUDE_DIRS})
  link_directories(${HDF5_LIBRARY_DIRS})
  set(LIBRARIES
    ${LIBRARIES}
    ${HDF5_LIBRARIES}
    )
endif(HDF5_FOUND)

find_package(MPI REQUIRED)

find_package(Qhull REQUIRED)
if (QHULL_FOUND)
  include_directories(${QHULL_INCLUDE_DIRS})
  link_directories(${QHULL_LIBRARY_DIRS})
  MESSAGE("-- QHULL include directory: ${QHULL_INCLUDE_DIRS}")
  MESSAGE("-- QHULL library directory: ${QHULL_LIBRARY}")
  MESSAGE("-- QHULL libraries: ${QHULL_LIBRARIES}")
  set(LIBRARIES
    ${LIBRARIES}
    ${QHULL_LIBRARIES}
    )
endif (QHULL_FOUND)

if (NOT QHULL_FOUND)
   MESSAGE(FATAL_ERROR "QHull not found")
endif ()

if (HAVE_OpenMP)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
endif(HAVE_OpenMP)

#   tools

option (PERF "detailed time measurement" OFF)
if (PERF)
  add_definitions("-DPERF")
endif(PERF)

#   fancy stuff

option (COLOROUTPUT "use colored output to the terminal" OFF)
if (COLOROUTPUT)
  add_definitions("-DCOLOROUTPUT")
endif(COLOROUTPUT)

option (DSERROR_DUMP "dserror creates a core file" OFF)
if (DSERROR_DUMP)
  add_definitions("-DDSERROR_DUMP")
endif(DSERROR_DUMP)

option (TRAP_FE "crash baci if a nan of inf occurs" ON)
if (TRAP_FE)
  add_definitions("-DTRAP_FE")
endif(TRAP_FE)

option (THROWELEMENTERRORS "throw errors at element level so all element error can be found" OFF)
if (THROWELEMENTERRORS)
  add_definitions("-DTHROWELEMENTERRORS")
endif(THROWELEMENTERRORS)

# Run some system checks
FILE(GLOB _check_files "cmake/checks/*.cmake")
FOREACH(_file ${_check_files})
  MESSAGE(STATUS "Include ${_file}")
  INCLUDE(${_file})
ENDFOREACH()

###------------------------------------------------------------------ External Libraries

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

#
### TRILINOS
SET(CMAKE_PREFIX_PATH ${Trilinos_PREFIX} ${CMAKE_PREFIX_PATH})
find_package(Trilinos REQUIRED)
# old find_package string for Trilinos < Q1/2012
#find_package(Trilinos REQUIRED CONFIGS TrilinosConfig_install.cmake)  # needed for build system on navier??

# Echo trilinos build info just for fun
MESSAGE("\nFound Trilinos!  Here are the details: ")
MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("End of Trilinos details\n")

# test for git executable
FIND_PROGRAM(Git_EXECUTABLE git
  DOC "git command line client")
MARK_AS_ADVANCED(Git_EXECUTABLE)

IF(Git_EXECUTABLE)
  EXEC_PROGRAM(${Git_EXECUTABLE} ${Trilinos_DIR} ARGS "describe" OUTPUT_VARIABLE Trilinos_GIT_VERSION )
  #MESSAGE(${Trilinos_GIT_VERSION})
  string (REGEX REPLACE "fatal.*$" "" Trilinos_GIT_VERSION "${Trilinos_GIT_VERSION}")
  ADD_DEFINITIONS( -DTrilinos_SHA1="${Trilinos_GIT_VERSION}")
ENDIF(Git_EXECUTABLE)

set(BOOST_LIBRARYDIR,${BOOST_LIBRARYDIR})
set(BOOST_INCLUDEDIR,${BOOST_INCLUDEDIR})

# This is needed on kaiser to enforce the usage of the non-system boost libraries
option (USE_NO_SYSTEM_PATHS "Avoid boost libraries in system paths" OFF)
#message("USE_NO_SYSTEM_PATHS: ${USE_NO_SYSTEM_PATHS}")
if (USE_NO_SYSTEM_PATHS)
  set(Boost_NO_SYSTEM_PATHS ON )
endif(USE_NO_SYSTEM_PATHS)

FIND_PACKAGE(Boost COMPONENTS graph program_options system filesystem REQUIRED)

if (Boost_FOUND)
  message("-- Found Boost: ${Boost_INCLUDE_DIRS}")
  message("-- Boost component libraries: ${Boost_LIBRARIES}")
  message("-- Boost libraries directory: ${Boost_LIBRARY_DIRS}")
  include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
  link_directories(${Boost_LIBRARY_DIRS})
  set(LIBRARIES
    ${Boost_LIBRARIES}
    ${LIBRARIES}
    )
  # Pass the version to baci
  add_definitions("-DBOOST_MAJOR_VERSION=${Boost_MAJOR_VERSION}")
  add_definitions("-DBOOST_MINOR_VERSION=${Boost_MINOR_VERSION}")

  # check Boost_LIBRARY_DIRS against Boost_INCLUDE_DIRS
  GET_FILENAME_COMPONENT(MYLIBPATH ${Boost_LIBRARY_DIRS} PATH)
  GET_FILENAME_COMPONENT(MYINCLUDEPATH ${Boost_INCLUDE_DIRS} PATH)
  string(COMPARE EQUAL "${MYLIBPATH}" "${MYINCLUDEPATH}" result)
  if(NOT result)
    message("WW ***********************************************************************")
    message("WW Warning: Paths for Boost library and include headers are inconsistent!")
    message("WW Warning: Please consult the CMakeLists.txt maintainer.")
    message("WW ***********************************************************************")
  endif()

  CHECK_CXX_SOURCE_COMPILES(
  "
  #include <boost/static_assert.hpp>
  int main()
  {
    BOOST_STATIC_ASSERT(2 < 3) __attribute__((unused));
  }
  "
  BOOST_STATIC_ASSERT_UNUSED_WORKS
  )

  if (BOOST_STATIC_ASSERT_UNUSED_WORKS)
    add_definitions(-DBACI_ATTRIBUTE_UNUSED=__attribute__\(\(unused\)\))
  else ()
    add_definitions(-DBACI_ATTRIBUTE_UNUSED=)
  endif()

else()
 message(FATAL_ERROR ${Boost_ERROR_REASON})
endif()

# trilinos libraries
# Choose trilinos libraries explicitly.
if (Trilinos_FOUND)
  REMOVE_DUPLICATES(Trilinos_LIBRARIES)
  REMOVE_DUPLICATES(Trilinos_TPL_LIBRARIES)
  REMOVE_DUPLICATES(Trilinos_INCLUDE_DIRS)
  REMOVE_DUPLICATES(Trilinos_TPL_INCLUDE_DIRS)

  # convert list of trilinos libraries to full path for consistency with other libraries
  SET(_libraries "")
  FOREACH(_library ${Trilinos_LIBRARIES})
    FIND_LIBRARY(TRILINOS_LIBRARY_${_library}
      NAMES ${_library}
      HINTS ${Trilinos_LIBRARY_DIRS}
      NO_DEFAULT_PATH
      NO_CMAKE_ENVIRONMENT_PATH
      NO_CMAKE_PATH
      NO_SYSTEM_ENVIRONMENT_PATH
      NO_CMAKE_SYSTEM_PATH
      NO_CMAKE_FIND_ROOT_PATH
      )
    LIST(APPEND _libraries ${TRILINOS_LIBRARY_${_library}})
  ENDFOREACH()

  include_directories(SYSTEM ${Trilinos_INCLUDE_DIRS})
  include_directories(SYSTEM ${Trilinos_TPL_INCLUDE_DIRS})
  link_directories(${Trilinos_LIBRARY_DIRS})

  set(LIBRARIES
    ${_libraries}
    ${Trilinos_TPL_LIBRARIES}
    ${LIBRARIES}
    )
endif (Trilinos_FOUND)

find_package(Netcdf REQUIRED)
if (NETCDF_FOUND)
  include_directories(${NETCDF_INCLUDE_DIRS})
  link_directories(${NETCDF_LIBRARY_DIRS})
  set(LIBRARIES
    ${LIBRARIES}
    ${NETCDF_LIBRARIES}
    )
endif(NETCDF_FOUND)

# deal.II package for deal.II support
CMAKE_POLICY(SET CMP0015 OLD)
find_package(deal.II 8.3 QUIET
    HINTS ${DEAL_II_DIR} ${deal.II_DIR} $ENV{DEAL_II_DIR}
  )
if (deal.II_FOUND)
  MESSAGE("\nFound deal.II!  Here are the details: ")
  MESSAGE("   deal.II_VERSION = ${DEAL_II_VERSION}")
  MESSAGE("   deal.II_DIR = ${DEAL_II_DIR}")
  include_directories(${DEAL_II_INCLUDE_DIRS})
  if (${CMAKE_BUILD_TYPE} STREQUAL "DEBUG")
    MESSAGE("   deal.II_LIBRARIES = ${DEAL_II_LIBRARIES_DEBUG}")
    link_directories(${DEAL_II_LIBRARIES_DEBUG})
    set(LIBRARIES
      ${LIBRARIES}
      ${DEAL_II_LIBRARIES_DEBUG}
      )
  else()
    MESSAGE("   deal.II_LIBRARIES = ${DEAL_II_LIBRARIES_RELEASE}")
    link_directories(${DEAL_II_LIBRARIES_RELEASE})
    set(LIBRARIES
      ${LIBRARIES}
      ${DEAL_II_LIBRARIES_RELEASE}
      )
  endif()
  MESSAGE("End of deal.II details\n")
  add_definitions("-DHAVE_DEAL_II")   # HAVE_DEAL_II preprocessor flag for BACI
else()
  message("WW *************************************************************")
  message("WW Warning: DEAL.II not found. No DEAL.II support in BACI.")
  message("WW Warning: You cannot run simulations that need DEAL.II.")
  message("WW *************************************************************")
endif()

# FFTW
find_package(FFTW)
if (FFTW_FOUND)
  MESSAGE("-- FFTW include directory: ${FFTW_INCLUDE_DIRS}")
  MESSAGE("-- FFTW library directory: ${FFTW_LIBRARY}")
  MESSAGE("-- FFTW libraries: ${FFTW_LIBRARIES}")
  include_directories(${FFTW_INCLUDE_DIRS})
  link_directories(${FFTW_LIBRARY_DIRS})
  set(LIBRARIES ${LIBRARIES} ${FFTW_LIBRARIES} )
  add_definitions("-DHAVE_FFTW")  # HAVE_FFTW preprocessor flag for BACI
  set(HAVE_FFTW ON)
else()
  message("WW *************************************************************")
  message("WW Warning: FFTW not found. no FFTW support in BACI.")
  message("WW Warning: You cannot run simulations that need a discrete FFT.")
  message("WW *************************************************************")
endif()

find_package(ACML)
if (ACML_FOUND)
  include_directories(${ACML_INCLUDES})
  #link_directories(${ACML_INCLUDES})
  set(LIBRARIES ${ACML_LIBRARIES} ${LIBRARIES}  )
  #set(LIBRARIES /lnm/lib/Q3_2012/amdlibm-3-0-2/lib/static/libamdlibm.a ${LIBRARIES} )
  MESSAGE("II Found ACML (AMD Core Math Library): ${ACML_LIBRARIES}")
  MESSAGE("II Found ACML (AMD Core Math Library): ${ACML_INCLUDES}")
endif()

# find AMD replacement for libM
find_package(AMDLIBM)
if (AMDLIBM_FOUND)
  include_directories(${AMDLIBM_INCLUDES})
  set(LIBRARIES ${AMDLIBM_LIBRARIES} ${LIBRARIES}  )
  MESSAGE("II Found AMD libM: ${AMDLIBM_LIBRARIES}")
  MESSAGE("II Found AMD libM: ${AMDLIBM_INCLUDES}")
endif()


find_package(CLN REQUIRED)
if (CLN_FOUND)
  include_directories(SYSTEM ${CLN_INCLUDE_DIRS})
  set(LIBRARIES ${CLN_LIBRARIES} ${LIBRARIES}  )
  MESSAGE("-- CLN include directory: ${CLN_INCLUDE_DIRS}")
  MESSAGE("-- CLN library directory: ${CLN_LIBRARY}")
  MESSAGE("-- CLN libraries: ${CLN_LIBRARIES}")
endif()
#explicitly remove all duplicate libraries, keeping the rightmost library
#additionally remove the keywords debug;optimized;general since the libraries
#are precompiled and this flag has no influence
REMOVE_DUPLICATES(LIBRARIES)
list(REMOVE_ITEM LIBRARIES debug optimized general)

if (HAVE_ASCII_FILE_COMPARISON_TESTS)
  message("-- Test comparing restults to ASCII files: enabled")
else()
  message("-- Test comparing restults to ASCII files: disabled")
endif()

if (HAVE_ENSIGHT_POSTPROCESSING_TESTS)
  message("-- Post-processing test for ensight filter: enabled")
else()
  message("-- Post-processing test for ensight filter: disabled")
endif()

###------------------------------------------------------------------ List of files
set(OBJS_DRT_ALE_LIB
  src/drt_ale2/ale2.cpp
  src/drt_ale2/ale2_evaluate.cpp
  src/drt_ale2/ale2_input.cpp
  src/drt_ale2/ale2_line.cpp
  src/drt_ale2/ale2_line_evaluate.cpp
  src/drt_ale2/ale2_nurbs.cpp
  src/drt_ale3/ale3.cpp
  src/drt_ale3/ale3_evaluate.cpp
  src/drt_ale3/ale3_input.cpp
  src/drt_ale3/ale3_nurbs.cpp
  src/drt_ale3/ale3_surface.cpp
  src/drt_ale3/ale3_surface_evaluate.cpp
)
set(OBJS_DRT_LIB
  src/drt_lib/drt_periodicbc.cpp
  src/drt_lib/drt_assemblestrategy.cpp
  src/drt_lib/drt_condition.cpp
  src/drt_lib/drt_conditiondefinition.cpp
  src/drt_lib/drt_container.cpp
  src/drt_lib/drt_discret.cpp
  src/drt_lib/drt_discret_faces.cpp
  src/drt_lib/drt_discret_hdg.cpp
  src/drt_lib/drt_discret_xwall.cpp
  src/drt_lib/drt_discret_xfem.cpp
  src/drt_lib/drt_discret_conditions.cpp
  src/drt_lib/drt_discret_evaluate.cpp
  src/drt_lib/drt_discret_fillcomplete.cpp
  src/drt_lib/drt_discret_partition.cpp
  src/drt_lib/drt_discret_utils.cpp
  src/drt_lib/drt_dofset.cpp
  src/drt_lib/drt_dofset_base.cpp
  src/drt_lib/drt_dofset_proxy.cpp
  src/drt_lib/drt_dofset_gidbased_wrapper.cpp
  src/drt_lib/drt_dofset_merged_wrapper.cpp
  src/drt_lib/drt_dofset_definedmapping_wrapper.cpp
  src/drt_lib/drt_dofset_transparent.cpp
  src/drt_lib/drt_dofset_transparent_independent.cpp
  src/drt_lib/drt_dofset_pbc.cpp
  src/drt_lib/drt_dofset_independent.cpp
  src/drt_lib/drt_dofset_independent_pbc.cpp
  src/drt_lib/drt_element.cpp
  src/drt_lib/drt_singletondestruction.cpp
  src/drt_lib/drt_function.cpp
  src/drt_lib/drt_functionvariables.cpp
  src/drt_lib/drt_globalproblem.cpp
  src/drt_lib/drt_locsys.cpp
  src/drt_lib/drt_node.cpp
  src/drt_lib/drt_parobject.cpp
  src/drt_lib/drt_parobjectfactory.cpp
  src/drt_lib/drt_utils_createdis.cpp
  src/drt_lib/drt_utils_materials.cpp
  src/drt_lib/drt_utils_factory.cpp
  src/drt_lib/drt_utils_nullspace.cpp
  src/drt_lib/drt_utils_discret.cpp
  src/drt_lib/drt_utils_discret_dbc.cpp
  src/drt_lib/drt_utils_discret_evaluate.cpp
  src/drt_lib/connectivity.cpp
)
set(OBJS_DRT_MAT_LIB
  src/drt_matelast/elast_aniso_structuraltensor_strategy.cpp
  src/drt_mat/aaa_mixedeffects.cpp
  src/drt_mat/aaagasser.cpp
  src/drt_mat/aaaneohooke.cpp
  src/drt_mat/aaaneohooke_stopro.cpp
  src/drt_mat/aaaraghavanvorp_damage.cpp
  src/drt_mat/arrhenius_pv.cpp
  src/drt_mat/arrhenius_spec.cpp
  src/drt_mat/arrhenius_temp.cpp
  src/drt_mat/beam_elasthyper.cpp
  src/drt_mat/beam_elasthyper_parameter.cpp
  src/drt_mat/carreauyasuda.cpp
  src/drt_mat/cavitationfluid.cpp
  src/drt_mat/cnst_1d_art.cpp
  src/drt_mat/compogden.cpp
  src/drt_mat/consolidation.cpp
  src/drt_mat/constraintmixture.cpp
  src/drt_mat/crosslinkermat.cpp
  src/drt_mat/constraintmixture_history.cpp
  src/drt_mat/damage.cpp
  src/drt_mat/elasthyper.cpp
  src/drt_mat/elasthyper_service.cpp
  src/drt_mat/elchmat.cpp
  src/drt_mat/elchphase.cpp
  src/drt_mat/elchsinglemat.cpp
  src/drt_mat/electrode.cpp
  src/drt_mat/electromagnetic.cpp
  src/drt_mat/ferech_pv.cpp
  src/drt_mat/fluid_linear_density_viscosity.cpp 
  src/drt_mat/fluid_murnaghantait.cpp
  src/drt_mat/fluidporo.cpp
  src/drt_mat/fluidporo_singlephase.cpp
  src/drt_mat/fluidporo_singlephaseDof.cpp
  src/drt_mat/fluidporo_singlephaselaw.cpp
  src/drt_mat/fluidporo_multiphase.cpp
  src/drt_mat/fluidporo_multiphase_reactions.cpp
  src/drt_mat/fluidporo_multiphase_singlereaction.cpp
  src/drt_mat/fluidporo_relpermeability_law.cpp
  src/drt_mat/fluidporo_viscosity_law.cpp
  src/drt_mat/fourieriso.cpp
  src/drt_mat/fouriervar.cpp
  src/drt_mat/membrane_elasthyper.cpp
  src/drt_mat/membrane_active_strain.cpp
  src/drt_mat/growthremodel_elasthyper.cpp
  src/drt_mat/growth.cpp
  src/drt_mat/growth_law.cpp
  src/drt_mat/herschelbulkley.cpp
  src/drt_mat/ion.cpp
  src/drt_mat/inelastic_defgrad_factors.cpp
  src/drt_mat/lubrication_law.cpp
  src/drt_mat/lubrication_mat.cpp
  src/drt_mat/material.cpp
  src/drt_mat/material_service.cpp
  src/drt_mat/matlist.cpp
  src/drt_mat/matlist_reactions.cpp
  src/drt_mat/matlist_bondreacs.cpp
  src/drt_mat/matpar_bundle.cpp
  src/drt_mat/matpar_material.cpp
  src/drt_mat/matpar_parameter.cpp
  src/drt_mat/micromaterial.cpp
  src/drt_mat/mixfrac.cpp
  src/drt_mat/modpowerlaw.cpp
  src/drt_mat/multiplicative_split_defgrad_elasthyper.cpp
  src/drt_mat/myocard.cpp
  src/drt_mat/myocard_fitzhugh_nagumo.cpp
  src/drt_mat/myocard_inada.cpp
  src/drt_mat/myocard_minimal.cpp
  src/drt_mat/myocard_san_garny.cpp
  src/drt_mat/myocard_tentusscher.cpp
  src/drt_mat/myocard_tools.cpp
  src/drt_mat/neohooke.cpp
  src/drt_mat/newman.cpp
  src/drt_mat/newman_multiscale.cpp
  src/drt_mat/newtonianfluid.cpp
  src/drt_mat/optimization_density.cpp
  src/drt_mat/permeablefluid.cpp
  src/drt_mat/plasticelasthyper.cpp
  src/drt_mat/plastic_VarConstUpdate.cpp
  src/drt_mat/plasticlinelast.cpp
  src/drt_mat/plasticnlnlogneohooke.cpp
  src/drt_mat/superelastic_sma.cpp
  src/drt_mat/poro_law.cpp
  src/drt_mat/poro_density_law.cpp
  src/drt_mat/robinson.cpp
  src/drt_mat/scalardepinterp.cpp
  src/drt_mat/scatra_mat.cpp
  src/drt_mat/scatra_mat_aniso.cpp
  src/drt_mat/scatra_mat_multiscale.cpp
  src/drt_mat/scatra_mat_poro_ecm.cpp
  src/drt_mat/scatra_mat_multiporo.cpp
  src/drt_mat/scatra_mat_var_chemdiffusion.cpp
  src/drt_mat/scatra_multiscale.cpp
  src/drt_mat/scatra_reaction_mat.cpp
  src/drt_mat/scatra_bondreac_mat.cpp
  src/drt_mat/scatra_reaction_coupling.cpp
  src/drt_mat/scatra_chemotaxis_mat.cpp
  src/drt_mat/matlist_reactions.cpp
  src/drt_mat/matlist_chemotaxis.cpp
  src/drt_mat/matlist_chemoreac.cpp
  src/drt_mat/scatra_mat_aniso.cpp
  src/drt_mat/soret.cpp
  src/drt_mat/spring.cpp
  src/drt_mat/structporo.cpp
  src/drt_mat/structporo_reaction.cpp
  src/drt_mat/structporo_reaction_ecm.cpp
  src/drt_mat/stvenantkirchhoff.cpp
  src/drt_mat/sutherland.cpp
  src/drt_mat/tempdepwater.cpp
  src/drt_mat/thermoplastichyperelast.cpp
  src/drt_mat/thermoplasticlinelast.cpp
  src/drt_mat/thermostvenantkirchhoff.cpp
  src/drt_mat/viscoanisotropic.cpp
  src/drt_mat/viscoelasthyper.cpp
  src/drt_mat/visconeohooke.cpp
  src/drt_mat/yoghurt.cpp
  src/drt_matelast/elast_anisoactivestress_evolution.cpp
  src/drt_matelast/elast_coupanisoexpoactive.cpp
  src/drt_matelast/elast_coupanisoexpo.cpp
  src/drt_matelast/elast_coupanisoexpotwocoup.cpp
  src/drt_matelast/elast_coupanisoneohooke.cpp
  src/drt_matelast/elast_coupanisoneohooke_VarProp.cpp
  src/drt_matelast/elast_coupanisopow.cpp
  src/drt_matelast/elast_couptransverselyisotropic.cpp
  src/drt_matelast/elast_coupblatzko.cpp
  src/drt_matelast/elast_couplogneohooke.cpp
  src/drt_matelast/elast_coupSaintVenantKirchhoff.cpp
  src/drt_matelast/elast_coupsimopister.cpp
  src/drt_matelast/elast_couplogmixneohooke.cpp
  src/drt_matelast/elast_coupexppol.cpp
  src/drt_matelast/elast_coupneohooke.cpp
  src/drt_matelast/elast_coupvarga.cpp
  src/drt_matelast/elast_isoanisoexpo.cpp
  src/drt_matelast/elast_isoexpopow.cpp
  src/drt_matelast/elast_coupmooneyrivlin.cpp
  src/drt_matelast/elast_isomooneyrivlin.cpp
  src/drt_matelast/elast_isoneohooke.cpp
  src/drt_matelast/elast_iso1pow.cpp
  src/drt_matelast/elast_iso2pow.cpp
  src/drt_matelast/elast_coup1pow.cpp
  src/drt_matelast/elast_coup2pow.cpp
  src/drt_matelast/elast_coup3pow.cpp
  src/drt_matelast/elast_coup13apow.cpp
  src/drt_matelast/elast_isovarga.cpp
  src/drt_matelast/elast_isoyeoh.cpp
  src/drt_matelast/elast_remodelfiber.cpp
  src/drt_matelast/elast_summand.cpp
  src/drt_matelast/elast_vologden.cpp
  src/drt_matelast/elast_volpenalty.cpp
  src/drt_matelast/elast_volsussmanbathe.cpp
  src/drt_matelast/elast_volpow.cpp
  src/drt_matelast/visco_isoratedep.cpp
  src/drt_matelast/visco_coupmyocard.cpp
  src/drt_matelast/visco_genmax.cpp
  src/drt_matelast/visco_fract.cpp
  src/drt_matelast/elast_isovolaaagasser.cpp
  src/drt_matelast/elast_isovolHUdependentneohooke.cpp
  src/drt_matelast/elast_isotestmaterial.cpp
  src/drt_patspec/patspec.cpp
  src/drt_mat/maxwell_0d_acinus.cpp
  src/drt_mat/maxwell_0d_acinus_NeoHookean.cpp
  src/drt_mat/maxwell_0d_acinus_Exponential.cpp
  src/drt_mat/maxwell_0d_acinus_DoubleExponential.cpp
  src/drt_mat/maxwell_0d_acinus_Ogden.cpp
  src/drt_mat/hemoglobin_0d_O2_saturation.cpp
  src/drt_mat/air_0d_O2_saturation.cpp
  src/drt_mat/particle_mat.cpp
  src/drt_mat/particle_mat_ellipsoids.cpp
  src/drt_mat/particle_material_base.cpp
  src/drt_mat/particle_material_thermo.cpp
  src/drt_mat/particle_material_sph_fluid.cpp
  src/drt_mat/particle_material_sph_boundary.cpp
  src/drt_mat/particle_material_dem.cpp
  src/drt_mat/acoustic.cpp
  src/drt_mat/acoustic_sol.cpp
  src/drt_mat/activefiber.cpp
  src/drt_mat/biochemo_mechano_cell_activefiber.cpp
  src/drt_mat/biochemo_mechano_cell_passivefiber.cpp
  src/drt_matelast/visco_generalizedgenmax.cpp
  src/drt_mat/mixture_elasthyper.cpp)

set(OBJS_DRT_ADAPT
  src/drt_adapter/adapter_algorithmbase.cpp
  src/drt_adapter/adapter_coupling.cpp
  src/drt_adapter/adapter_coupling_mortar.cpp
  src/drt_adapter/adapter_coupling_nonlin_mortar.cpp
  src/drt_adapter/adapter_coupling_poro_mortar.cpp
  src/drt_adapter/adapter_coupling_ehl_mortar.cpp
  src/drt_adapter/adapter_coupling_volmortar.cpp
  src/drt_adapter/ad_field_wrapper.cpp
  src/drt_adapter/ad_fld_fluid.cpp
  src/drt_adapter/ad_fld_fluid_ale.cpp
  src/drt_adapter/ad_fld_fluid_ale_immersed.cpp
  src/drt_adapter/ad_fld_base_algorithm.cpp
  src/drt_adapter/ad_fld_lung.cpp
  src/drt_adapter/ad_fld_moving_boundary.cpp
  src/drt_adapter/ad_fld_poro.cpp
  src/drt_adapter/ad_fld_fluid_fsi_msht.cpp
  src/drt_adapter/ad_fld_fluid_fsi.cpp
  src/drt_adapter/ad_fld_fluid_ac_fsi.cpp
  src/drt_adapter/ad_fld_fluid_fpsi.cpp
  src/drt_adapter/ad_fld_fluid_fluid_fsi.cpp
  src/drt_adapter/ad_fld_fluid_xfem.cpp
  src/drt_adapter/ad_fld_fluid_immersed.cpp
  src/drt_adapter/ad_fld_fluid_ale_xfem.cpp
  src/drt_adapter/ad_fld_fluid_xfsi.cpp
  src/drt_adapter/ad_fld_fluid_fluid_xfsi.cpp
  src/drt_adapter/ad_opt.cpp
  src/drt_adapter/ad_opt_fluid_adjoint_base.cpp
  src/drt_adapter/ad_opt_fluid_algo.cpp
  src/drt_adapter/adapter_lubrication.cpp
  src/drt_adapter/adapter_scatra_base_algorithm.cpp
  src/drt_adapter/adapter_scatra_fluid_ale_coupling_algo.cpp
  src/drt_adapter/adapter_scatra_fluid_coupling_algorithm.cpp
  src/drt_adapter/adapter_particle.cpp
  src/drt_adapter/ad_str_structure.cpp
  src/drt_structure/stru_dyn_nln_drt.cpp
  src/drt_adapter/ad_str_factory.cpp
  src/drt_adapter/ad_str_structure_new.cpp
  src/drt_adapter/ad_str_loca.cpp
  src/drt_adapter/ad_str_loca_wrapper.cpp
  src/drt_adapter/ad_str_constr_merged.cpp
  src/drt_adapter/ad_str_lung.cpp
  src/drt_adapter/ad_str_redairway.cpp
  src/drt_adapter/ad_str_timint_adaptive.cpp
  src/drt_adapter/ad_str_wrapper.cpp
  src/drt_adapter/ad_str_ssiwrapper.cpp
  src/drt_adapter/ad_str_fsiwrapper.cpp
  src/drt_adapter/ad_str_fpsiwrapper.cpp
  src/drt_adapter/ad_str_pasiwrapper.cpp
  src/drt_adapter/ad_str_structalewrapper.cpp
  src/drt_adapter/ad_str_fsiwrapper_immersed.cpp
  src/drt_adapter/ad_str_fsi_timint_adaptive.cpp
  src/drt_adapter/ad_str_poro_wrapper.cpp
  src/drt_adapter/ad_str_timeloop.cpp
  src/drt_adapter/ad_str_xcontact.cpp
  src/drt_adapter/ad_str_xstructure_algorithm.cpp
  src/drt_adapter/ad_str_invana.cpp
  src/drt_adapter/adapter_thermo.cpp
  src/drt_adapter/ad_ale.cpp
  src/drt_adapter/ad_ale_fluid.cpp
  src/drt_adapter/ad_ale_fpsi.cpp
  src/drt_adapter/ad_ale_fsi_msht.cpp
  src/drt_adapter/ad_ale_fsi.cpp
  src/drt_adapter/ad_ale_wear.cpp
  src/drt_adapter/ad_ale_wrapper.cpp
  src/drt_adapter/ad_ale_xffsi.cpp
  src/drt_adapter/ad_porofluidmultiphase_wrapper.cpp
  src/drt_adapter/ad_str_multiphysicswrapper_cellmigration.cpp
  src/drt_adapter/ad_scatra_wrapper.cpp
  src/drt_adapter/ad_scatra_wrapper_cellmigration.cpp
  
)
set(OBJS_DRT_CONTACT_LIB
  src/drt_contact/contact_analytical.cpp
  src/drt_contact/contact_element.cpp
  src/drt_contact/contact_node.cpp
  src/drt_contact/friction_node.cpp
)
set(OBJS_DRT_RED_AIRWAYS_LIB
  src/drt_red_airways/airway_evaluate.cpp
  src/drt_red_airways/acinus_evaluate.cpp
  src/drt_red_airways/inter_acinar_dep_evaluate.cpp
  src/drt_red_airways/red_air_blood_scatra_evaluate.cpp
  src/drt_red_airways/red_air_blood_scatraLine3_evaluate.cpp
  src/drt_red_airways/airway_impl.cpp
  src/drt_red_airways/acinus_impl.cpp
  src/drt_red_airways/inter_acinar_dep_impl.cpp
  src/drt_red_airways/red_air_blood_scatra_impl.cpp
  src/drt_red_airways/red_air_blood_scatraLine3_impl.cpp
  src/drt_red_airways/red_airway.cpp
  src/drt_red_airways/red_acinus.cpp
  src/drt_red_airways/red_inter_acinar_dep.cpp
  src/drt_red_airways/red_air_blood_scatra.cpp
  src/drt_red_airways/red_air_blood_scatraLine3.cpp
  src/drt_red_airways/red_airway_input.cpp
  src/drt_red_airways/red_airway_resulttest.cpp
)
set(OBJS_DRT_INPAR
  src/drt_inpar/drt_boolifyparameters.cpp
  src/drt_inpar/drt_validconditions.cpp
  src/drt_inpar/drt_validmaterials.cpp
  src/drt_inpar/drt_validparameters.cpp
  src/drt_inpar/inpar_cardiovascular0d.cpp
  src/drt_inpar/inpar_acou.cpp
  src/drt_inpar/inpar_ale.cpp
  src/drt_inpar/inpar_beamcontact.cpp
  src/drt_inpar/inpar_beaminteraction.cpp
  src/drt_inpar/inpar_beampotential.cpp
  src/drt_inpar/inpar_browniandyn.cpp
  src/drt_inpar/inpar_geometry_pair.cpp
  src/drt_inpar/inpar_bio.cpp
  src/drt_inpar/inpar_cardiac_monodomain.cpp
  src/drt_inpar/inpar_particle.cpp
  src/drt_inpar/inpar_particle_old.cpp
  src/drt_inpar/inpar_cavitation.cpp
  src/drt_inpar/inpar_pasi.cpp
  src/drt_inpar/inpar_cell.cpp
  src/drt_inpar/inpar_contact.cpp
  src/drt_inpar/inpar_cut.cpp
  src/drt_inpar/inpar_ehl.cpp
  src/drt_inpar/inpar_elemag.cpp
  src/drt_inpar/inpar_elch.cpp
  src/drt_inpar/inpar_fluid.cpp
  src/drt_inpar/inpar_fpsi.cpp
  src/drt_inpar/inpar_fs3i.cpp
  src/drt_inpar/inpar_fsi.cpp
  src/drt_inpar/inpar_immersed.cpp
  src/drt_inpar/inpar_invanalysis.cpp
  src/drt_inpar/inpar_io.cpp
  src/drt_inpar/inpar_IO_monitor_structure_dbc.cpp
  src/drt_inpar/inpar_IO_runtime_vtk_output.cpp
  src/drt_inpar/inpar_IO_runtime_vtk_output_structure.cpp
  src/drt_inpar/inpar_IO_runtime_vtk_output_structure_beams.cpp
  src/drt_inpar/inpar_IO_runtime_vtp_output_structure.cpp
  src/drt_inpar/inpar_levelset.cpp
  src/drt_inpar/inpar_loca_continuation.cpp
  src/drt_inpar/inpar_lubrication.cpp
  src/drt_inpar/inpar_binningstrategy.cpp
  src/drt_inpar/inpar_mlmc.cpp
  src/drt_inpar/inpar_mor.cpp
  src/drt_inpar/inpar_mortar.cpp
  src/drt_inpar/inpar_poroelast.cpp
  src/drt_inpar/inpar_plasticity.cpp
  src/drt_inpar/inpar_poroscatra.cpp
  src/drt_inpar/inpar_poromultiphase.cpp
  src/drt_inpar/inpar_poromultiphase_scatra.cpp
  src/drt_inpar/inpar_porofluidmultiphase.cpp
  src/drt_inpar/inpar_problemtype.cpp
  src/drt_inpar/inpar_s2i.cpp
  src/drt_inpar/inpar_scatra.cpp
  src/drt_inpar/inpar_searchtree.cpp
  src/drt_inpar/inpar_solver.cpp
  src/drt_inpar/inpar_solver_nonlin.cpp
  src/drt_inpar/inpar_ssi.cpp
  src/drt_inpar/inpar_sti.cpp
  src/drt_inpar/inpar_structure.cpp
  src/drt_inpar/inpar_thermo.cpp
  src/drt_inpar/inpar_topopt.cpp
  src/drt_inpar/inpar_tsi.cpp
  src/drt_inpar/inpar_twophase.cpp
  src/drt_inpar/inpar_volmortar.cpp
  src/drt_inpar/inpar_wear.cpp
  src/drt_inpar/inpar_contact_xcontact.cpp
  src/drt_inpar/inpar_xfem.cpp
  src/drt_tutorial/inpar_tutorial.cpp
)
set(OBJS_DRT_ARTNET
  src/drt_art_net/art_net_dyn_drt.cpp
  src/drt_art_net/art_net_utils.cpp
  src/drt_art_net/art_net_timint.cpp
  src/drt_art_net/art_net_impl_stationary.cpp
  src/drt_art_net/artnetexplicitintegration.cpp
)


set(OBJS_DRT_XFEM
  src/drt_xfem/xfem_edgestab.cpp
  src/drt_xfem/xfem_discretization_utils.cpp
  src/drt_xfem/xfem_dofset_transparent_independent.cpp
  src/drt_xfem/xfem_dofset.cpp
  src/drt_xfem/xfem_mesh_projector.cpp
  src/drt_xfem/xfem_neumann.cpp
  src/drt_xfem/xfem_interface_utils.cpp
  src/drt_xfem/xfem_utils.cpp
  src/drt_xfem/xfem_state_creator.cpp
  src/drt_xfem/xfield_state.cpp
  src/drt_xfem/xfield_field_coupling_dofset.cpp
  src/drt_xfem/xfield_field_coupling.cpp
  src/drt_xfem/xfem_multi_field_mapextractor.cpp
  src/drt_xfem/xfluid_timeInt_base.cpp
  src/drt_xfem/xfluid_timeInt_std_SemiLagrange.cpp
  src/drt_xfem/xfluid_timeInt.cpp
)


set(OBJS_DRT_XFEM_MANAG_LIB
  src/drt_xfem/xfem_coupling_levelset.cpp
  src/drt_xfem/xfem_coupling_base.cpp
  src/drt_xfem/xfem_coupling_mesh.cpp
  src/drt_xfem/xfem_coupling_fpi_mesh.cpp
  src/drt_xfem/xfem_coupling_mesh_coupled_levelset.cpp
  src/drt_xfem/xfem_condition_manager.cpp
)

set(OBJS_POST_DRT_MONITOR
  src/post_monitor/post_drt_monitor.cpp
)
set(OBJS_DRT_ARTNET_LIB
  src/drt_art_net/artery.cpp
  src/drt_art_net/artery_ele_calc_lin_exp.cpp
  src/drt_art_net/artery_ele_calc_pres_based.cpp
  src/drt_art_net/artery_ele_factory.cpp
  src/drt_art_net/artery_evaluate.cpp
  src/drt_art_net/artery_ele_calc.cpp
  src/drt_art_net/art_junction.cpp
  src/drt_art_net/art_terminal_bc.cpp
  src/drt_art_net/art_write_gnuplot.cpp
  src/drt_art_net/artery_input.cpp
  src/drt_art_net/artery_resulttest.cpp
)
set(OBJS_FILTER_COMMON
  src/filter_common/filter_evaluation.cpp
)
set(OBJS_DRT_BEAMCONTACT_LIB
  src/drt_beamcontact/beam3contact_manager.cpp
  src/drt_beamcontact/beam3contact_octtree.cpp
  src/drt_beamcontact/beam3contactnew.cpp
  src/drt_beamcontact/beam3contactinterface.cpp
  src/drt_beamcontact/beam3contact.cpp
  src/drt_beamcontact/beam3contactvariables.cpp
  src/drt_beamcontact/beam3tospherecontact.cpp
  src/drt_beamcontact/beam3tosolidcontact.cpp
  src/drt_beamcontact/beam3tosolidmeshtying.cpp
  src/drt_beamcontact/beam3tobeampotential.cpp
  src/drt_beamcontact/beam3tospherepotential.cpp
  src/drt_beamcontact/str_model_evaluator_beaminteraction_old.cpp
)
set(OBJS_DRT_BEAMINTERACTION_LIB
  src/drt_beaminteraction/periodic_boundingbox.cpp
  src/drt_beaminteraction/beam3contact_utils.cpp
  src/drt_beaminteraction/beam3contact_tangentsmoothing.cpp
  src/drt_beaminteraction/beam_contact_params.cpp
  src/drt_beaminteraction/beam_contact_runtime_vtk_output_params.cpp
  src/drt_beaminteraction/beam_contact_pair.cpp
  src/drt_beaminteraction/beam_to_beam_contact_pair.cpp
  src/drt_beaminteraction/beam_to_beam_contact_params.cpp
  src/drt_beaminteraction/beam_to_beam_contact_variables.cpp
  src/drt_beaminteraction/beam_to_solid_mortar_manager.cpp
  src/drt_beaminteraction/beam_to_solid_volume_meshtying_pair_base.cpp
  src/drt_beaminteraction/beam_to_solid_volume_meshtying_pair_gauss_point.cpp
  src/drt_beaminteraction/beam_to_solid_volume_meshtying_pair_gauss_point_cross_section.cpp
  src/drt_beaminteraction/beam_to_solid_volume_meshtying_pair_mortar.cpp
  src/drt_beaminteraction/beam_to_solid_volume_meshtying_pair_factory.cpp
  src/drt_beaminteraction/beam_to_solid_volume_meshtying_params.cpp
  src/drt_beaminteraction/beam_to_solid_volume_meshtying_vtk_output_params.cpp
  src/drt_beaminteraction/beam_to_solid_volume_meshtying_vtk_output_writer.cpp
  src/drt_beaminteraction/beam_to_solid_vtu_output_writer_base.cpp
  src/drt_beaminteraction/beam_to_solid_vtu_output_writer_visualization.cpp
  src/drt_beaminteraction/beam_to_sphere_contact_pair.cpp
  src/drt_beaminteraction/beam_to_sphere_contact_params.cpp
  src/drt_beaminteraction/beam_potential_params.cpp
  src/drt_beaminteraction/beam_potential_runtime_vtk_output_params.cpp
  src/drt_beaminteraction/beam_potential_pair.cpp
  src/drt_beaminteraction/beam_to_beam_potential_pair.cpp
  src/drt_beaminteraction/beam_to_sphere_potential_pair.cpp
  src/drt_beaminteraction/beam_link.cpp
  src/drt_beaminteraction/beam_link_rigidjointed.cpp
  src/drt_beaminteraction/beam_link_beam3r_lin2_rigidjointed.cpp
  src/drt_beaminteraction/beam_link_pinjointed.cpp
  src/drt_beaminteraction/beam_link_beam3r_lin2_pinjointed.cpp
  src/drt_beaminteraction/beam_link_truss.cpp
  src/drt_beaminteraction/str_model_evaluator_beamdata.cpp
  src/drt_beaminteraction/str_model_evaluator_beaminteraction.cpp
  src/drt_beaminteraction/str_model_evaluator_beaminteraction_datastate.cpp
  src/drt_beaminteraction/beaminteraction_submodel_evaluator_factory.cpp
  src/drt_beaminteraction/beaminteraction_submodel_evaluator_generic.cpp
  src/drt_beaminteraction/beaminteraction_submodel_evaluator_beamcontact.cpp
  src/drt_beaminteraction/beaminteraction_submodel_evaluator_beamcontact_assembly_manager.cpp
  src/drt_beaminteraction/beaminteraction_submodel_evaluator_beamcontact_assembly_manager_direct.cpp
  src/drt_beaminteraction/beaminteraction_submodel_evaluator_beamcontact_assembly_manager_indirect.cpp
  src/drt_beaminteraction/beaminteraction_submodel_evaluator_potential.cpp
  src/drt_beaminteraction/beaminteraction_calc_utils.cpp
  src/drt_beaminteraction/beaminteraction_data.cpp
  src/drt_beaminteraction/crosslinking_params.cpp
  src/drt_beaminteraction/beaminteraction_submodel_evaluator_crosslinking.cpp
  src/drt_beaminteraction/spherebeamlinking_params.cpp
  src/drt_beaminteraction/beaminteraction_submodel_evaluator_spherebeamlinking.cpp
  src/drt_beaminteraction/beaminteraction_geometry_utils.cpp
)
set(OBJS_DRT_STRU_MULTI
  src/drt_stru_multi/convert_mat.cpp
  src/drt_stru_multi/microstatic.cpp
  src/drt_stru_multi/microstatic_npsupport.cpp
  src/drt_stru_multi/microstatic_service.cpp
)
set(OBJS_SOH8_LIB
  src/drt_so3/so_utils.cpp
  src/drt_so3/inversedesign.cpp
  src/drt_so3/inversedesign_evaluate.cpp
  src/drt_so3/prestress.cpp
  src/drt_so3/so_base.cpp
  src/drt_so3/so_hex18.cpp
  src/drt_so3/so_hex20.cpp
  src/drt_so3/so_hex20_evaluate.cpp
  src/drt_so3/so_hex20_input.cpp
  src/drt_so3/so_hex20_multiscale.cpp
  src/drt_so3/so_hex20_service.cpp
  src/drt_so3/so_hex27.cpp
  src/drt_so3/so_hex27_evaluate.cpp
  src/drt_so3/so_hex27_input.cpp
  src/drt_so3/so_hex27_multiscale.cpp
  src/drt_so3/so_hex27_service.cpp
  src/drt_so3/so_hex8.cpp
  src/drt_so3/so_hex8_eas.cpp
  src/drt_so3/so_hex8_evaluate.cpp
  src/drt_so3/so_hex8_input.cpp
  src/drt_so3/so_hex8_multiscale.cpp
  src/drt_so3/so_hex8_service.cpp
  src/drt_so3/so_hex8_determinant_analysis.cpp
  src/drt_so3/so_hex8p1j1.cpp
  src/drt_so3/so_hex8p1j1_evaluate.cpp
  src/drt_so3/so_hex8p1j1_input.cpp
  src/drt_so3/so_hex8fbar.cpp
  src/drt_so3/so_hex8fbar_evaluate.cpp
  src/drt_so3/so_hex8fbar_input.cpp
  src/drt_so3/so_line.cpp
  src/drt_so3/so_line_evaluate.cpp
  src/drt_so3/so_nstet5.cpp
  src/drt_so3/so_nstet5_evaluate.cpp
  src/drt_so3/so_nstet5_input.cpp
  src/drt_so3/so_nstet5_nodalstrain.cpp
  src/drt_so3/so_nstet5_multiscale.cpp
  src/drt_so3/so_nstet.cpp
  src/drt_so3/so_nstet_evaluate.cpp
  src/drt_so3/so_nstet_input.cpp
  src/drt_so3/so_nstet_nodalstrain.cpp
  src/drt_so3/so_nstet_multiscale.cpp
  src/drt_so3/so_nurbs27.cpp
  src/drt_so3/so_nurbs27_evaluate.cpp
  src/drt_so3/so_nurbs27_input.cpp
  src/drt_so3/so_sh18.cpp
  src/drt_so3/so_sh18_input.cpp
  src/drt_so3/so_sh18_evaluate.cpp
  src/drt_so3/so_sh8.cpp
  src/drt_so3/so_sh8_evaluate.cpp
  src/drt_so3/so_sh8_input.cpp
  src/drt_so3/so_sh8_service.cpp
  src/drt_so3/so_sh8p8.cpp
  src/drt_so3/so_sh8p8_auxiliar.cpp
  src/drt_so3/so_sh8p8_eas.cpp
  src/drt_so3/so_sh8p8_evaluate.cpp
  src/drt_so3/so_sh8p8_input.cpp
  src/drt_so3/so_shw6.cpp
  src/drt_so3/so_shw6_evaluate.cpp
  src/drt_so3/so_shw6_input.cpp
  src/drt_so3/so_surface.cpp
  src/drt_so3/so_surface_evaluate.cpp
  src/drt_so3/so_surface_traceestimate.cpp
  src/drt_so3/so_tet10.cpp
  src/drt_so3/so_tet10_evaluate.cpp
  src/drt_so3/so_tet10_input.cpp
  src/drt_so3/so_tet4.cpp
  src/drt_so3/so_tet4_evaluate.cpp
  src/drt_so3/so_tet4_multiscale.cpp
  src/drt_so3/so_tet4_input.cpp
  src/drt_so3/so_tet4av.cpp
  src/drt_so3/so_tet4av_evaluate.cpp
  src/drt_so3/so_tet4av_input.cpp
  src/drt_so3/so_weg6.cpp
  src/drt_so3/so_weg6_evaluate.cpp
  src/drt_so3/so_weg6_input.cpp
  src/drt_so3/so_pyramid5.cpp
  src/drt_so3/so_pyramid5_evaluate.cpp
  src/drt_so3/so_pyramid5_input.cpp
  src/drt_so3/so_pyramid5_service.cpp
  src/drt_so3/so_pyramid5_multiscale.cpp
  src/drt_so3/so_pyramid5fbar.cpp
  src/drt_so3/so_pyramid5fbar_evaluate.cpp
  src/drt_so3/so_pyramid5fbar_input.cpp
  src/drt_so3/so3_poro.cpp
  src/drt_so3/so3_poro_evaluate.cpp
  src/drt_so3/so3_poro_eletypes.cpp
  src/drt_so3/so3_poro_p1.cpp
  src/drt_so3/so3_poro_p1_eletypes.cpp
  src/drt_so3/so3_poro_p1_evaluate.cpp
  src/drt_so3/so3_poro_p1_scatra.cpp
  src/drt_so3/so3_poro_p1_scatra_eletypes.cpp
  src/drt_so3/so3_poro_scatra.cpp
  src/drt_so3/so3_poro_scatra_eletypes.cpp
  src/drt_so3/so3_scatra.cpp
  src/drt_so3/so3_scatra_evaluate.cpp
  src/drt_so3/so3_scatra_eletypes.cpp
  src/drt_so3/so3_thermo.cpp
  src/drt_so3/so3_thermo_evaluate.cpp
  src/drt_so3/so3_thermo_eletypes.cpp
  src/drt_so3/so3_plast/so3_ssn_plast.cpp
  src/drt_so3/so3_plast/so3_ssn_plast_evaluate.cpp
  src/drt_so3/so3_plast/so3_ssn_plast_eletypes.cpp
  src/drt_so3/so3_plast/so3_ssn_plast_eas.cpp
  src/drt_so3/so3_plast/so3_ssn_plast_sosh8.cpp
  src/drt_so3/so3_plast/so3_ssn_plast_sosh18.cpp
)
set(OBJS_DRT_ELCH
  src/drt_elch/elch_algorithm.cpp
  src/drt_elch/elch_dyn.cpp
  src/drt_elch/elch_moving_boundary_algorithm.cpp
)
set(OBJS_DRT_VAR
  src/drt_variational/var_chemdiff_dyn.cpp
)
set(OBJS_DRT_BELE3_LIB
  src/drt_bele3/bele2.cpp
  src/drt_bele3/bele2_evaluate.cpp
  src/drt_bele3/bele3.cpp
  src/drt_bele3/bele3_evaluate.cpp
  src/drt_bele3/bele3_line.cpp
  src/drt_bele3/bele3_line_evaluate.cpp
  src/drt_bele3/vele3.cpp
  src/drt_bele3/vele3_evaluate.cpp
  src/drt_bele3/vele3_line.cpp
  src/drt_bele3/vele3_line_evaluate.cpp
  src/drt_bele3/vele3_surface.cpp
  src/drt_bele3/vele3_surface_evaluate.cpp
)
set(OBJS_DRT_COMM
  src/drt_comm/comm_utils.cpp
)

set(OBJS_DRT_UTILS
  src/drt_binstrategy/binning_strategy_utils.cpp
  src/drt_binstrategy/binning_strategy.cpp
  src/drt_binstrategy/drt_meshfree_bin.cpp
  src/drt_binstrategy/drt_meshfree_multibin.cpp
  src/drt_fem_general/debug_nan.cpp
  src/drt_fem_general/drt_utils_boundary_integration.cpp
  src/drt_fem_general/drt_utils_fem_shapefunctions.cpp
  src/drt_fem_general/drt_utils_integration.cpp
  src/drt_fem_general/drt_utils_gausspoints.cpp
  src/drt_fem_general/drt_utils_local_connectivity_matrices.cpp
  src/drt_fem_general/drt_utils_polynomial.cpp
  src/drt_fem_general/drt_utils_shapevalues_hdg.cpp
  src/drt_fem_general/largerotations.cpp
  src/drt_fluid/fluid_functions.cpp
  src/drt_fluid/fluid_rotsym_periodicbc_utils.cpp
  src/drt_fluid_xfluid/xfluid_functions.cpp
  src/drt_fluid_xfluid/xfluid_functions_combust.cpp
  src/drt_poromultiphase_scatra/poromultiphase_scatra_function.cpp
  src/drt_geometry/integrationcell.cpp
  src/drt_geometry/intersection_interfacepoint.cpp
  src/drt_geometry/intersection_math.cpp
  src/drt_geometry/intersection_service.cpp
  src/drt_geometry/position_array.cpp
  src/drt_geometry/searchtree.cpp
  src/drt_geometry/searchtree_geometry_service.cpp
  src/drt_geometry/searchtree_nearestobject.cpp
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
  src/drt_lib/drt_utils_metis.cpp
  src/drt_lib/drt_utils_parallel.cpp
  src/drt_lib/epetra_utils.cpp
  src/drt_beaminteraction/drt_utils_parallel_proctoproc.cpp
  src/drt_lib/drt_utils_parmetis.cpp
  src/pss_full/pss_am.c
  src/pss_full/pss_table.c
  src/pss_full/pss_table_iter.c
)
set(OBJS_DRT_ERROR
  src/drt_lib/drt_dserror.cpp
)
set(OBJS_DRT_EHL
  src/drt_ehl/ehl_base.cpp
  src/drt_ehl/ehl_dyn.cpp
  src/drt_ehl/ehl_monolithic.cpp
  src/drt_ehl/ehl_partitioned.cpp
  src/drt_ehl/ehl_utils.cpp
)
set(OBJS_DRT_LUBRICATION
  src/drt_lubrication/lubrication_dyn.cpp
  src/drt_lubrication/lubrication_resulttest.cpp
  src/drt_lubrication/lubrication_timint_implicit.cpp
  src/drt_lubrication/lubrication_timint_stat.cpp
)
set(OBJS_DRT_SCATRA
  src/drt_scatra/scatra_algorithm.cpp
  src/drt_scatra/scatra_dyn.cpp
  src/drt_scatra/scatra_cardiac_monodomain_dyn.cpp
  src/drt_levelset/levelset_algorithm.cpp
  src/drt_levelset/levelset_algorithm_utils.cpp
  src/drt_levelset/levelset_algorithm_reinit.cpp
  src/drt_levelset/levelset_dyn.cpp
  src/drt_scatra/scatra_resulttest.cpp
  src/drt_scatra/scatra_resulttest_elch.cpp
  src/drt_scatra/scatra_resulttest_hdg.cpp
  src/drt_scatra/scatra_timint_bdf2.cpp
  src/drt_scatra/scatra_timint_cardiac_monodomain_scheme.cpp
  src/drt_scatra/scatra_timint_cardiac_monodomain.cpp
  src/drt_scatra/scatra_timint_convcheck_strategies.cpp
  src/drt_scatra/scatra_timint_genalpha.cpp
  src/drt_scatra/scatra_timint_implicit.cpp
  src/drt_scatra/scatra_timint_implicit_service.cpp
  src/drt_scatra/scatra_timint_ost.cpp
  src/drt_scatra/scatra_timint_stat.cpp
  src/drt_scatra/scatra_timint_loma.cpp
  src/drt_scatra/scatra_timint_loma_genalpha.cpp
  src/drt_scatra/scatra_timint_loma_ost.cpp
  #src/drt_scatra/scatra_timint_loma_bdf2.cpp
  src/drt_scatra/scatra_timint_meshtying_strategy_fluid.cpp
  src/drt_scatra/scatra_timint_heterogeneous_reaction_strategy.cpp
  src/drt_scatra/scatra_timint_meshtying_strategy_fluid_elch.cpp
  src/drt_scatra/scatra_timint_meshtying_strategy_s2i.cpp
  src/drt_scatra/scatra_timint_meshtying_strategy_s2i_elch.cpp
  src/drt_scatra/scatra_timint_meshtying_strategy_std.cpp
  src/drt_scatra/scatra_timint_meshtying_strategy_std_elch.cpp
  src/drt_scatra/scatra_timint_meshtying_strategy_std_var_chemdiff.cpp 
  src/drt_scatra/scatra_timint_meshtying_strategy_artery.cpp
  src/drt_scatra/scatra_timint_elch.cpp
  src/drt_scatra/scatra_timint_elch_scheme.cpp
  src/drt_scatra/scatra_timint_poromulti.cpp
  src/drt_scatra/scatra_timint_variational.cpp
  src/drt_scatra/scatra_timint_variational_scheme.cpp 
  src/drt_scatra/scatra_utils_clonestrategy.cpp
  src/drt_scatra/scatra_utils.cpp
  src/drt_scatra/turbulence_hit_initial_scalar_field.cpp
  src/drt_scatra/turbulence_hit_scalar_forcing.cpp
  src/drt_levelset/levelset_intersection_utils.cpp
  src/drt_levelset/levelset_timint_ost.cpp
  src/drt_levelset/levelset_timint_stat.cpp
  src/drt_scatra/scatra_timint_hdg.cpp
  src/drt_scatra/scatra_timint_cardiac_monodomain_scheme_hdg.cpp
  src/drt_contact_xcontact/xcontact_levelset_timint_ost.cpp
)
set(OBJS_DRT_TSI
  src/drt_tsi/tsi_algorithm.cpp
  src/drt_tsi/tsi_monolithic.cpp
  src/drt_tsi/tsi_partitioned.cpp
  src/drt_tsi/tsi_dyn.cpp
  src/drt_tsi/tsi_utils.cpp
)
set(OBJS_DRT_S8_LIB
  src/drt_s8/shell8.cpp
  src/drt_s8/shell8_evaluate.cpp
  src/drt_s8/shell8_input.cpp
  src/drt_s8/shell8_line.cpp
  src/drt_s8/shell8_line_evaluate.cpp
  src/drt_s8/shell8_scatra.cpp
  src/fortran/s8_fserv.f
  src/drt_s8/shell8/s8_mat_linel.c
  src/drt_s8/shell8/s8_mat_ogden.c
  src/drt_s8/shell8/s8_mat_ogden2.c
  src/drt_s8/shell8/s8_mat_ogden3.c
  src/drt_s8/shell8/s8_mattransform.c
  src/drt_s8/shell8/s8_tmat.c
)
set(OBJS_DRT_MEMBRANE_LIB
  src/drt_membrane/membrane_eletypes.cpp
  src/drt_membrane/membrane.cpp
  src/drt_membrane/membrane_evaluate.cpp
  src/drt_membrane/membrane_input.cpp
  src/drt_membrane/membrane_line.cpp
  src/drt_membrane/membrane_line_evaluate.cpp
  src/drt_membrane/membrane_scatra_eletypes.cpp
  src/drt_membrane/membrane_scatra.cpp
  src/drt_membrane/membrane_scatra_evaluate.cpp
)
set(OBJS_DRT_FSI
  src/drt_fsi/fsi_str_model_evaluator_partitioned.cpp
  src/drt_fsi/fs_monolithic.cpp
  src/drt_fsi/fsi_algorithm.cpp
  src/drt_fsi/fsi_constrmonolithic.cpp
  src/drt_fsi/fsi_constrmonolithic_fluidsplit.cpp
  src/drt_fsi/fsi_constrmonolithic_structuresplit.cpp
  src/drt_fsi/fsi_debugwriter.cpp
  src/drt_fsi/fsi_dirichletneumann_factory.cpp
  src/drt_fsi/fsi_dirichletneumann.cpp
  src/drt_fsi/fsi_dirichletneumann_disp.cpp
  src/drt_fsi/fsi_dirichletneumannslideale.cpp
  src/drt_fsi/fsi_dirichletneumann_volcoupl.cpp
  src/drt_fsi/fsi_dyn.cpp
  src/drt_fsi/fsi_fluid_ale.cpp
  src/drt_fsi/fsi_lungmonolithic.cpp
  src/drt_fsi/fsi_lungmonolithic_fluidsplit.cpp
  src/drt_fsi/fsi_lungmonolithic_structuresplit.cpp
  src/drt_fsi/fsi_matrixtransform.cpp
  src/drt_fsi/fsi_monolithic.cpp
  src/drt_fsi/fsi_monolithic_redistribution.cpp
  src/drt_fsi/fsi_monolithic_timintada.cpp
  src/drt_fsi/fsi_monolithic_nonox.cpp
  src/drt_fsi/fsi_monolithic_linearsystem.cpp
  src/drt_fsi/fsi_monolithicfluidsplit.cpp
  src/drt_fsi/fsi_monolithicstructuresplit.cpp
  src/drt_fsi/fsi_fluidfluidmonolithic_structuresplit_nonox.cpp
  src/drt_fsi/fsi_fluidfluidmonolithic_structuresplit.cpp
  src/drt_fsi/fsi_fluidfluidmonolithic_fluidsplit_nonox.cpp
  src/drt_fsi/fsi_fluidfluidmonolithic_fluidsplit.cpp
  src/drt_fsi/fsi_mortarmonolithic_structuresplit.cpp
  src/drt_fsi/fsi_mortarmonolithic_fluidsplit.cpp
  src/drt_fsi/fsi_nox_aitken.cpp
  src/drt_fsi/fsi_nox_epsilon.cpp
  src/drt_fsi/fsi_nox_fixpoint.cpp
  src/drt_fsi/fsi_nox_group.cpp
  src/drt_fsi/fsi_nox_jacobian.cpp
  src/drt_fsi/fsi_nox_linearsystem_gcr.cpp
  src/drt_fsi/fsi_nox_mpe.cpp
  src/drt_fsi/fsi_nox_newton.cpp
  src/drt_fsi/fsi_nox_sd.cpp
  src/drt_fsi/fsi_overlapprec.cpp
  src/drt_fsi/fsi_overlapprec_fsiamg.cpp
  src/drt_fsi/fsi_overlapprec_hybrid.cpp
  src/drt_fsi/fsi_lung_overlapprec.cpp
  src/drt_fsi/fsi_constr_overlapprec.cpp
  src/drt_fsi/fsi_overlapprec_fsiamg_analysis.cpp
  src/drt_fsi/fsi_overlapprec_fsiamg_smoothers.cpp
  src/drt_fsi/fsi_overlapprec_fsiamg_vcycles.cpp
  src/drt_fsi/fsi_partitioned.cpp
  src/drt_fsi/fsi_resulttest.cpp
  src/drt_fsi/fsi_slidingmonolithic_fluidsplit.cpp
  src/drt_fsi/fsi_slidingmonolithic_structuresplit.cpp
  src/drt_fsi/fsi_statustest.cpp
  src/drt_fsi/fsi_structureale.cpp
  src/drt_fsi/fsi_utils.cpp
  src/drt_fsi/fsi_overlapprec_amgnxn.cpp
)
set(OBJS_DRT_FSI_XFEM
  src/drt_fsi_xfem/fsi_xfem_fluid.cpp
  src/drt_fsi_xfem/fsi_xfem_monolithic.cpp
  src/drt_fsi_xfem/fsi_xfem_algorithm.cpp
  src/drt_fsi_xfem/coupling_comm_manager.cpp
  src/drt_fsi_xfem/XFScoupling_manager.cpp
  src/drt_fsi_xfem/XFAcoupling_manager.cpp
  src/drt_fsi_xfem/XFPcoupling_manager.cpp
  src/drt_fsi_xfem/XFFcoupling_manager.cpp
)
set(OBJS_DRT_FS3I
  src/drt_fs3i/fs3i.cpp
  src/drt_fs3i/fs3i_dyn.cpp
  src/drt_fs3i/fs3i_partitioned.cpp
  src/drt_fs3i/fps3i_partitioned.cpp
  src/drt_fs3i/fs3i_partitioned_1wc.cpp
  src/drt_fs3i/fps3i_partitioned_1wc.cpp
  src/drt_fs3i/fs3i_partitioned_2wc.cpp
  src/drt_fs3i/ac_fsi.cpp
  src/drt_fs3i/ac_fsi_lts.cpp
  src/drt_fs3i/biofilm_fsi.cpp
  src/drt_fs3i/biofilm_fsi_utils.cpp
)
set(OBJS_DRT_RED_AIRWAYS
  src/drt_red_airways/airwayimplicitintegration.cpp
  src/drt_red_airways/redairway_tissue.cpp
  src/drt_red_airways/red_airways_dyn_drt.cpp
)
set(OBJS_DRT_THERMO_LIB
  src/drt_thermo/thermo_ele_boundary_impl.cpp
  src/drt_thermo/thermo_ele_impl.cpp
  src/drt_thermo/thermo_ele_impl_utils.cpp
  src/drt_thermo/thermo_element.cpp
  src/drt_thermo/thermo_element_boundary_evaluate.cpp
  src/drt_thermo/thermo_element_evaluate.cpp
  src/drt_thermo/thermo_element_input.cpp
)
set(OBJS_DRT_REGISTER
  src/drt_lib/drt_parobjectregister.cpp
)
set(OBJS_DRT_LUBRICATION_LIB
  src/drt_lubrication_ele/lubrication_ele.cpp
  src/drt_lubrication_ele/lubrication_ele_evaluate.cpp
  src/drt_lubrication_ele/lubrication_ele_factory.cpp
  src/drt_lubrication_ele/lubrication_ele_calc.cpp
  src/drt_lubrication_ele/lubrication_ele_parameter.cpp

)
set(OBJS_DRT_POROFLUIDMULTIPHASE_LIB
  src/drt_porofluidmultiphase_ele/porofluidmultiphase_ele.cpp
  src/drt_porofluidmultiphase_ele/porofluidmultiphase_ele_calc.cpp
  src/drt_porofluidmultiphase_ele/porofluidmultiphase_ele_calc_utils.cpp
  src/drt_porofluidmultiphase_ele/porofluidmultiphase_ele_factory.cpp
  src/drt_porofluidmultiphase_ele/porofluidmultiphase_ele_parameter.cpp
  src/drt_porofluidmultiphase_ele/porofluidmultiphase_ele_evaluate.cpp
  src/drt_porofluidmultiphase_ele/porofluid_phasemanager.cpp
  src/drt_porofluidmultiphase_ele/porofluid_variablemanager.cpp
  src/drt_porofluidmultiphase_ele/porofluid_evaluator.cpp
  src/drt_porofluidmultiphase_ele/porofluidmultiphase_ele_boundary_evaluate.cpp
  src/drt_porofluidmultiphase_ele/porofluidmultiphase_ele_boundary_factory.cpp
  src/drt_porofluidmultiphase_ele/porofluidmultiphase_ele_boundary_calc.cpp
)
set(OBJS_DRT_SCATRA_LIB
  src/drt_scatra_ele/scatra_ele.cpp
  src/drt_scatra_ele/scatra_ele_boundary_calc.cpp
  src/drt_scatra_ele/scatra_ele_boundary_calc_elch.cpp
  src/drt_scatra_ele/scatra_ele_boundary_calc_elch_diffcond.cpp
  src/drt_scatra_ele/scatra_ele_boundary_calc_elch_electrode.cpp
  src/drt_scatra_ele/scatra_ele_boundary_calc_elch_electrode_growth.cpp
  src/drt_scatra_ele/scatra_ele_boundary_calc_elch_electrode_sti_thermo.cpp
  src/drt_scatra_ele/scatra_ele_boundary_calc_elch_NP.cpp
  src/drt_scatra_ele/scatra_ele_boundary_calc_loma.cpp
  src/drt_scatra_ele/scatra_ele_boundary_calc_poro.cpp
  src/drt_scatra_ele/scatra_ele_boundary_calc_refconc_reac.cpp
  src/drt_scatra_ele/scatra_ele_boundary_calc_std.cpp
  src/drt_scatra_ele/scatra_ele_boundary_calc_sti_electrode.cpp
  src/drt_scatra_ele/scatra_ele_boundary_evaluate.cpp
  src/drt_scatra_ele/scatra_ele_boundary_factory.cpp
  src/drt_scatra_ele/scatra_ele_calc.cpp
  src/drt_scatra_ele/scatra_ele_calc_variational.cpp
  src/drt_scatra_ele/scatra_ele_calc_advanced_reaction.cpp
  src/drt_scatra_ele/scatra_ele_calc_refconc_reac.cpp
  src/drt_scatra_ele/scatra_ele_calc_aniso.cpp
  src/drt_scatra_ele/scatra_ele_calc_bondreac.cpp
  src/drt_scatra_ele/scatra_ele_calc_cardiac_monodomain.cpp
  src/drt_scatra_ele/scatra_ele_calc_chemo.cpp
  src/drt_scatra_ele/scatra_ele_calc_chemo_reac.cpp
  src/drt_scatra_ele/scatra_ele_calc_elch.cpp
  src/drt_scatra_ele/scatra_ele_calc_elch_diffcond.cpp
  src/drt_scatra_ele/scatra_ele_calc_elch_diffcond_multiscale.cpp
  src/drt_scatra_ele/scatra_ele_calc_elch_diffcond_sti_thermo.cpp
  src/drt_scatra_ele/scatra_ele_calc_elch_electrode.cpp
  src/drt_scatra_ele/scatra_ele_calc_elch_electrode_sti_thermo.cpp
  src/drt_scatra_ele/scatra_ele_calc_elch_NP.cpp
  src/drt_scatra_ele/scatra_ele_calc_loma.cpp
  src/drt_scatra_ele/scatra_ele_calc_ls.cpp
  src/drt_scatra_ele/scatra_ele_calc_lsreinit.cpp
  src/drt_scatra_ele/scatra_ele_calc_OD.cpp
  src/drt_scatra_ele/scatra_ele_calc_poro.cpp
  src/drt_scatra_ele/scatra_ele_calc_poro_reac.cpp
  src/drt_scatra_ele/scatra_ele_calc_multiporo_reac.cpp
  src/drt_scatra_ele/scatra_ele_calc_artery.cpp
  src/drt_scatra_ele/scatra_ele_calc_poro_reac_ECM.cpp
  src/drt_scatra_ele/scatra_ele_calc_service_general.cpp
  src/drt_scatra_ele/scatra_ele_calc_service_stabilization.cpp
  src/drt_scatra_ele/scatra_ele_calc_service_turbulence.cpp
  src/drt_scatra_ele/scatra_ele_calc_service_ls.cpp
  src/drt_scatra_ele/scatra_ele_calc_service_lsreinit.cpp
  src/drt_scatra_ele/scatra_ele_calc_service_loma.cpp
  src/drt_scatra_ele/scatra_ele_calc_service_elch_NP.cpp
  src/drt_scatra_ele/scatra_ele_calc_service_elch_electrode.cpp
  src/drt_scatra_ele/scatra_ele_calc_service_elch_diffcond.cpp
  src/drt_scatra_ele/scatra_ele_calc_service_elch_diffcond_multiscale.cpp
  src/drt_scatra_ele/scatra_ele_calc_service_elch.cpp
  src/drt_scatra_ele/scatra_ele_calc_service_cardiac_monodomain.cpp
  src/drt_scatra_ele/scatra_ele_calc_std.cpp
  src/drt_scatra_ele/scatra_ele_calc_sti_diffcond.cpp
  src/drt_scatra_ele/scatra_ele_calc_sti_electrode.cpp
  src/drt_scatra_ele/scatra_ele_calc_utils.cpp
  src/drt_scatra_ele/scatra_ele_evaluate.cpp
  src/drt_scatra_ele/scatra_ele_factory.cpp
  src/drt_scatra_ele/scatra_ele_input.cpp
  src/drt_scatra_ele/scatra_ele_parameter_elch.cpp
  src/drt_scatra_ele/scatra_ele_parameter_elch_diffcond.cpp
  src/drt_scatra_ele/scatra_ele_parameter_lsreinit.cpp
  src/drt_scatra_ele/scatra_ele_parameter_std.cpp
  src/drt_scatra_ele/scatra_ele_parameter_timint.cpp
  src/drt_scatra_ele/scatra_ele_parameter_turbulence.cpp
  src/drt_scatra_ele/scatra_ele_sti_elch.cpp
  src/drt_scatra_ele/scatra_ele_sti_thermo.cpp
  src/drt_scatra_ele/scatra_ele_utils_elch.cpp
  src/drt_scatra_ele/scatra_ele_utils_elch_diffcond.cpp
  src/drt_scatra_ele/scatra_ele_utils_elch_electrode.cpp
  src/drt_scatra_ele/scatra_ele_hdg.cpp
  src/drt_scatra_ele/scatra_ele_calc_hdg.cpp
  src/drt_scatra_ele/scatra_ele_hdg_intfaces_calc.cpp
  src/drt_scatra_ele/scatra_ele_hdg_boundary_calc.cpp
  src/drt_scatra_ele/scatra_ele_calc_cardiac_monodomain_hdg.cpp
  src/drt_scatra_ele/scatra_ele_calc_service_no_physics.cpp
  src/drt_scatra_ele/scatra_ele_calc_no_physics.cpp
)
set(OBJS_DRT_STI
  src/drt_sti/sti_algorithm.cpp
  src/drt_sti/sti_clonestrategy.cpp
  src/drt_sti/sti_dyn.cpp
  src/drt_sti/sti_monolithic.cpp
  src/drt_sti/sti_partitioned.cpp
  src/drt_sti/sti_resulttest.cpp
)
set(OBJS_DRT_THERMO
  src/drt_thermo/thr_aux.cpp
  src/drt_thermo/thr_dyn.cpp
  src/drt_thermo/thr_resulttest.cpp
  src/drt_thermo/thrtimint.cpp
  src/drt_thermo/thrtimint_genalpha.cpp
  src/drt_thermo/thrtimint_impl.cpp
  src/drt_thermo/thrtimint_ost.cpp
  src/drt_thermo/thrtimint_statics.cpp
  src/drt_thermo/thrtimint_expl.cpp
  src/drt_thermo/thrtimint_expleuler.cpp
)
set(OBJS_DRT_FLUID
  src/drt_fluid/fluid_coupling_red_models.cpp
  src/drt_fluid/fluid_dyn_nln_drt.cpp
  src/drt_fluid/fluid_meshtying.cpp
  src/drt_fluid/fluid_discret_extractor.cpp
  src/drt_fluid/fluid_utils.cpp
  src/drt_fluid/fluid_utils_infnormscaling.cpp
  src/drt_fluid/fluid_utils_mapextractor.cpp
  src/drt_fluid/fluid_MHD_evaluate.cpp
  src/drt_fluid/fluid_timint_bdf2.cpp
  src/drt_fluid/fluid_timint_genalpha.cpp
  src/drt_fluid/fluid_timint_loma_genalpha.cpp
  src/drt_fluid/fluid_timint_loma_ost.cpp
  src/drt_fluid/fluid_timint_loma_bdf2.cpp
  src/drt_fluid/fluid_timint_loma.cpp
  src/drt_fluid/fluid_timint_ost.cpp
  src/drt_fluid/fluid_timint_poro_genalpha.cpp
  src/drt_fluid/fluid_timint_poro_ost.cpp
  src/drt_fluid/fluid_timint_poro_stat.cpp
  src/drt_fluid/fluid_timint_poro.cpp
  src/drt_fluid/fluid_timint_ac_ost.cpp
  src/drt_fluid/fluid_timint_ac.cpp
  src/drt_fluid/fluid_timint_red_ost.cpp
  src/drt_fluid/fluid_timint_red_bdf2.cpp
  src/drt_fluid/fluid_timint_red_genalpha.cpp
  src/drt_fluid/fluid_timint_red_stat.cpp
  src/drt_fluid/fluid_timint_red.cpp
  src/drt_fluid/fluid_timint_stat.cpp
  src/drt_fluid/fluid_timint_topopt_bdf2.cpp
  src/drt_fluid/fluid_timint_topopt_genalpha.cpp
  src/drt_fluid/fluid_timint_topopt_ost.cpp
  src/drt_fluid/fluid_timint_topopt_stat.cpp
  src/drt_fluid/fluid_timint_topopt.cpp
  src/drt_fluid/fluid_timint_two_phase_genalpha.cpp
  src/drt_fluid/fluid_timint_two_phase_ost.cpp
  src/drt_fluid/fluid_timint_two_phase_stat.cpp
  src/drt_fluid/fluid_timint_two_phase.cpp
  src/drt_fluid/fluidimpedancecondition.cpp
  src/drt_fluid/fluid_volumetric_surfaceFlow_condition.cpp
  src/drt_fluid/fluidimplicitintegration.cpp
  src/drt_fluid/fluidresulttest.cpp
  src/drt_fluid/fluid_timint.cpp
  src/drt_fluid/fluid_timint_hdg.cpp
  src/drt_fluid/fluid_timint_stat_hdg.cpp
  src/drt_fluid/fluid_xwall.cpp
)
set(OBJS_DRT_FLUID_TURBULENCE
  src/drt_fluid_turbulence/dyn_smag.cpp
  src/drt_fluid_turbulence/dyn_vreman.cpp
  src/drt_fluid_turbulence/boxfilter.cpp
  src/drt_fluid_turbulence/drt_transfer_turb_inflow.cpp
  src/drt_fluid_turbulence/turbulence_hit_initial_field.cpp
  src/drt_fluid_turbulence/turbulence_hit_forcing.cpp
  src/drt_fluid_turbulence/turbulent_flow_algorithm.cpp
  src/drt_fluid_turbulence/turbulence_statistic_manager.cpp
  src/drt_fluid_turbulence/turbulence_statistics_bfs.cpp
  src/drt_fluid_turbulence/turbulence_statistics_ph.cpp
  src/drt_fluid_turbulence/turbulence_statistics_bfda.cpp
#  src/drt_fluid_turbulence/turbulence_statistics_bcf.cpp
  src/drt_fluid_turbulence/turbulence_statistics_ccy.cpp
  src/drt_fluid_turbulence/turbulence_statistics_cha.cpp
  src/drt_fluid_turbulence/turbulence_statistics_hit.cpp
  src/drt_fluid_turbulence/turbulence_statistics_tgv.cpp
  src/drt_fluid_turbulence/turbulence_statistics_ldc.cpp
  src/drt_fluid_turbulence/turbulence_statistics_mean_general.cpp
#  src/drt_fluid_turbulence/turbulence_statistics_oracles.cpp
  src/drt_fluid_turbulence/turbulence_statistics_sqc.cpp
)
set(OBJS_DRT_FLUID_XFLUID
  src/drt_fluid_xfluid/xfluid.cpp
  src/drt_fluid_xfluid/xfluidfluid.cpp
  src/drt_fluid_xfluid/xfluidfluid_state.cpp
  src/drt_fluid_xfluid/xfluidresulttest.cpp
  src/drt_fluid_xfluid/xfluid_state.cpp
  src/drt_fluid_xfluid/xfluid_state_creator.cpp
  src/drt_fluid_xfluid/xfluid_levelset_coupling_algorithm.cpp
  src/drt_fluid_xfluid/xfluid_outputservice.cpp
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
set(OBJS_DRT_CONSTRAINT_LIB
  src/drt_constraint/constraint_element2.cpp
  src/drt_constraint/constraint_element2_evaluate.cpp
  src/drt_constraint/constraint_element3.cpp
  src/drt_constraint/constraint_element3_evaluate.cpp
)
set(OBJS_DRT_NURBS_DISCRET_LIB
  src/drt_fem_general/drt_utils_bspline.cpp
  src/drt_fem_general/drt_utils_nurbs_shapefunctions.cpp
  src/drt_nurbs_discret/drt_control_point.cpp
  src/drt_nurbs_discret/drt_knotvector.cpp
  src/drt_nurbs_discret/drt_nurbs_discret.cpp
)
set(OBJS_DRT_BINSTRATEGY_DISCRET_LIB
  src/drt_beaminteraction/crosslinker_node.cpp
  src/drt_particle/particle_ellipsoid_node.cpp
  src/drt_particle/particle_node.cpp
  src/drt_particle/particle_radius_node.cpp
)
set(OBJS_DRT_TOPOPT
  src/drt_opti/topopt_algorithm.cpp
  src/drt_opti/topopt_dyn.cpp
  src/drt_opti/topopt_utils.cpp
)

set(OBJS_DRT_TOPOPT_LIB
  src/drt_opti/topopt_optimizer_ele.cpp
  src/drt_opti/topopt_optimizer_ele_evaluate.cpp
  src/drt_opti/topopt_optimizer_ele_impl.cpp
  src/drt_opti/topopt_optimizer_ele_parameter.cpp
)

set(OBJS_DRT_TOPOPT_ADJOINT
  src/drt_opti/topopt_fluidAdjointImplTimeIntegration.cpp
  src/drt_opti/topopt_fluidAdjointResulttest.cpp
  src/drt_opti/topopt_fluidAdjoint_timeint.cpp
)

set(OBJS_OPTI
  src/drt_opti/topopt_optimizer.cpp
  src/drt_opti/opti_GCMMA.cpp
  src/drt_opti/opti_resulttest.cpp
)

set(OBJS_DRT_LOMA
  src/drt_loma/loma_algorithm.cpp
  src/drt_loma/loma_dyn.cpp
)

set(OBJS_DRT_CAVITATION
  src/drt_particle/cavitation_algorithm.cpp
  src/drt_particle/cavitation_service.cpp
  src/drt_particle/particle_dyn.cpp
)

set(OBJS_DRT_TWO_PHASE
  src/drt_two_phase_flow/two_phase_algorithm.cpp
  src/drt_two_phase_flow/two_phase_dyn.cpp
)

set(OBJS_DRT_CONTACT_XCONTACT
  src/drt_contact_xcontact/xcontact_algorithm_base.cpp
  src/drt_contact_xcontact/xcontact_algorithm_partitioned.cpp
  src/drt_contact_xcontact/xcontact_dyn.cpp
  src/drt_contact_xcontact/xcontact_multi_discretization_wrapper.cpp
  src/drt_contact_xcontact/xcontact_levelset_algorithm.cpp
  src/drt_contact_xcontact/xcontact_levelset_intersection.cpp
  src/drt_contact_xcontact/xcontact_levelset_reinit_generic.cpp
  src/drt_contact_xcontact/xcontact_levelset_reinit_elliptic.cpp
  src/drt_contact_xcontact/xcontact_levelset_reinit_elliptic_nox_interface.cpp
)

set(OBJS_POST_DRT_COMMON
  src/post_drt_common/post_drt_common.cpp
)
set(OBJS_POST_PROCESSOR_MAIN
  src/post_drt_common/post_processor.cpp
  src/post_drt_common/post_single_field_writers.cpp
  src/post_drt_common/post_structure_stress.cpp
  src/post_drt_common/post_thermo_heatflux.cpp
  src/post_drt_common/post_filter_base.cpp
  src/post_drt_common/post_writer_base.cpp
)
set(OBJS_POST_DRT_ENSIGHT
  src/post_ensight/post_drt_ensight_writer.cpp
  src/post_ensight/post_drt_ensight_writer_nurbs.cpp
)
set(OBJS_POST_DRT_VTK
  src/post_vtk/post_drt_vtu_writer.cpp
  src/post_vtk/post_drt_vtu_writer_node_based.cpp
  src/post_vtk/post_drt_vti_writer.cpp
  src/post_vtk/post_drt_vtp_writer.cpp
  src/post_vtk/post_drt_vtk_writer.cpp
)
set(OBJS_INV_ANALYSIS
  src/drt_inv_analysis/gen_inv_analysis.cpp
  src/drt_structure/str_invanalysis.cpp
  src/drt_inv_analysis/invana_cal_drt.cpp
)
set(OBJS_STATINV_ANALYSIS
  src/drt_inv_analysis/invana_cal_drt.cpp
  src/drt_inv_analysis/invana_control.cpp
  src/drt_inv_analysis/invana_writer.cpp
  src/drt_inv_analysis/timint_adjoint.cpp
  src/drt_inv_analysis/timint_adjoint_prestress.cpp
  src/drt_inv_analysis/objective_funct_disp.cpp
  src/drt_inv_analysis/objective_funct_surfcurr.cpp
  src/drt_inv_analysis/matpar_manager.cpp
  src/drt_inv_analysis/matpar_manager_uniform.cpp
  src/drt_inv_analysis/matpar_manager_elementwise.cpp
  src/drt_inv_analysis/matpar_manager_patchwise.cpp
  src/drt_inv_analysis/matpar_manager_tvsvd.cpp
  src/drt_inv_analysis/regularization_base.cpp
  src/drt_inv_analysis/regularization_totalvariation.cpp
  src/drt_inv_analysis/regularization_tikhonov.cpp
  src/drt_inv_analysis/invana_resulttest.cpp
  src/drt_inv_analysis/invana_utils.cpp
  src/drt_inv_analysis/invana_factory.cpp
  src/drt_inv_analysis/invana_base.cpp
  src/drt_inv_analysis/invana_auglagr.cpp
  src/drt_inv_analysis/initial_guess.cpp
  src/drt_inv_analysis/optimizer_factory.cpp
  src/drt_inv_analysis/optimizer_base.cpp
  src/drt_inv_analysis/optimizer_lbfgs.cpp
  src/drt_inv_analysis/optimizer_smc.cpp
  src/drt_inv_analysis/optimizer_smc_predict.cpp
  src/drt_inv_analysis/optimizer_mh.cpp
  src/drt_inv_analysis/optimizer_bruteforce.cpp
  src/drt_inv_analysis/particle_group.cpp
  src/drt_inv_analysis/particle_data.cpp
  src/drt_inv_analysis/particle_comm.cpp
  src/drt_inv_analysis/likelihood_evaluation.cpp
  src/drt_inv_analysis/DcsMatrix.cpp
  src/drt_inv_analysis/chol_factory.cpp
  src/drt_inv_analysis/chol_factor.cpp
  src/drt_inv_analysis/metropolis_kernel.cpp
)
set(OBJS_STATINV_ANALYSIS_LIB
  src/drt_inv_analysis/particle_data.cpp
)
set(OBJS_MLMC
	src/drt_mlmc/drt_uq_redairways.cpp
	src/drt_mlmc/drt_uq_dyn.cpp
  src/drt_mlmc/mlmc.cpp
  src/drt_mlmc/mc_utils.cpp
  src/drt_mlmc/mc_mat_par_manager.cpp
  src/drt_mlmc/mc_var_thickness_manager.cpp
  src/drt_mlmc/randomfield.cpp
  src/drt_mlmc/randomfield_fourier.cpp
  src/drt_mlmc/randomfield_spectral.cpp
  src/drt_mlmc/randomvariable.cpp
)
set(OBJS_DRT_MAT
  src/drt_mat/micromaterial_evaluate.cpp
  src/drt_mat/micromaterialgp_static.cpp
)
set(OBJS_DRT_MAT_SCATRA_MULTI
  src/drt_mat/scatra_multiscale_evaluate.cpp
  src/drt_mat/scatra_multiscale_gp.cpp
)
set(OBJS_BACI_MAIN
  src/global_full/global_init_control.cpp
  src/global_full/global_inp_control.cpp
  src/global_full/global_cal_control.cpp
  src/global_full/global_control.cpp
  src/global_full/baci.cpp
)
set(OBJS_DRT_RIGIDSPHERE_LIB
  src/drt_rigidsphere/rigidsphere.cpp
  src/drt_rigidsphere/rigidsphere_evaluate.cpp
  src/drt_rigidsphere/rigidsphere_input.cpp
)
set(OBJS_PRE_LOCSYS
  src/pre_locsys/pre_locsys.cpp
)
set(OBJS_DRT_SURFSTRESS
  src/drt_surfstress/drt_surfstress_manager.cpp
)
set(OBJS_DRT_TO3_LIB
  src/drt_torsion3/torsion3.cpp
  src/drt_torsion3/torsion3_evaluate.cpp
  src/drt_torsion3/torsion3_input.cpp
)
set(OBJS_DRT_MORTAR
  src/drt_mortar/base_binarytree.cpp
  src/drt_mortar/mortar_binarytree.cpp
  src/drt_mortar/mortar_coupling2d.cpp
  src/drt_mortar/mortar_coupling3d.cpp
  src/drt_mortar/mortar_coupling3d_classes.cpp
  src/drt_mortar/mortar_utils.cpp
  src/drt_mortar/mortar_dofset.cpp
  src/drt_mortar/mortar_integrator.cpp
  src/drt_mortar/mortar_interface.cpp
  src/drt_mortar/mortar_interface_tools.cpp
  src/drt_mortar/mortar_manager_base.cpp
  src/drt_mortar/mortar_strategy_factory.cpp
  src/drt_mortar/mortar_projector.cpp
  src/drt_mortar/mortar_strategy_base.cpp
  src/drt_mortar/mortar_matrix_transform.cpp
  src/drt_mortar/mortar_multifield_coupling.cpp
)
set(OBJS_DRT_MORTAR_LIB
  src/drt_mortar/mortar_element.cpp
  src/drt_mortar/mortar_element_integrator.cpp
  src/drt_mortar/mortar_element_shapefct.cpp
  src/drt_mortar/mortar_node.cpp
  src/drt_contact/contact_nitsche_utils.cpp
)
set(OBJS_DRT_VOLMORTAR
src/drt_volmortar/volmortar_cell.cpp
src/drt_volmortar/volmortar_coupling.cpp
src/drt_volmortar/volmortar_integrator.cpp
src/drt_volmortar/volmortar_utils.cpp
)
set(OBJS_DRT_STRUCTURE
  src/drt_structure/strtimada.cpp
  src/drt_structure/strtimada_create.cpp
  src/drt_structure/strtimada_zienxie.cpp
  src/drt_structure/strtimint.cpp
  src/drt_structure/strtimint_create.cpp
  src/drt_structure/strtimint_ab2.cpp
  src/drt_structure/strtimint_centrdiff.cpp
  src/drt_structure/strtimint_expl.cpp
  src/drt_structure/strtimint_expleuler.cpp
  src/drt_structure/strtimint_gemm.cpp
  src/drt_structure/strtimint_genalpha.cpp
  src/drt_structure/strtimint_impl.cpp
  src/drt_structure/strtimint_impl_nox.cpp
  src/drt_structure/strtimint_noxgroup.cpp
  src/drt_structure/strtimint_noxlinsys.cpp
  src/drt_structure/strtimint_ost.cpp
  src/drt_structure/strtimint_statics.cpp
  src/drt_structure/strtimint_prestress.cpp
  src/drt_structure/stru_aux.cpp
  src/drt_structure/stru_resulttest.cpp
  src/drt_structure_new/str_utils.cpp
  src/drt_structure_new/str_dbc.cpp
  src/drt_structure_new/str_monitor_dbc.cpp
#  src/drt_beaminteraction/str_dbc_biopolynet.cpp
  src/drt_structure_new/str_resulttest.cpp
  src/drt_structure_new/str_factory.cpp
  src/drt_structure_new/str_impl_generic.cpp
  src/drt_structure_new/str_impl_statics.cpp
  src/drt_structure_new/str_impl_genalpha.cpp
  src/drt_structure_new/str_impl_genalpha_liegroup.cpp
  src/drt_structure_new/str_impl_gemm.cpp
  src/drt_structure_new/str_impl_ost.cpp
  src/drt_structure_new/str_integrator.cpp
  src/drt_structure_new/str_model_evaluator.cpp
  src/drt_structure_new/str_model_evaluator_data.cpp
  src/drt_structure_new/str_model_evaluator_contactdata.cpp
  src/drt_browniandyn/str_model_evaluator_browniandyn.cpp
  src/drt_browniandyn/str_model_evaluator_browniandyndata.cpp
  src/drt_structure_new/str_model_evaluator_factory.cpp
  src/drt_structure_new/str_model_evaluator_generic.cpp
  src/drt_structure_new/str_model_evaluator_structure.cpp
  src/drt_structure_new/str_model_evaluator_cardiovascular0d.cpp
  src/drt_structure_new/str_model_evaluator_lagpenconstraint.cpp
  src/drt_structure_new/str_model_evaluator_springdashpot.cpp
  src/drt_structure_new/str_model_evaluator_contact.cpp
  src/drt_structure_new/str_model_evaluator_meshtying.cpp
  src/drt_contact_xcontact/str_model_evaluator_xcontact.cpp
  src/drt_struct_ale/struct_ale_str_model_evaluator.cpp
  src/drt_structure_new/str_model_evaluator_multiphysics.cpp
  src/drt_immersed_problem/str_model_evaluator_multiphysics_cellmigration.cpp
  src/drt_immersed_problem/str_model_evaluator_partitioned_ssi_adhesiondynamics.cpp
  src/drt_structure_new/str_nln_solver_factory.cpp
  src/drt_structure_new/str_nln_solver_utils.cpp
  src/drt_structure_new/str_nln_solver_generic.cpp
  src/drt_structure_new/str_nln_solver_fullnewton.cpp
  src/drt_structure_new/str_nln_solver_nox.cpp
  src/drt_structure_new/str_nln_solver_ptc.cpp
#  src/drt_structure_new/str_nln_solver_nonlin.cpp
  src/drt_structure_new/str_predict_factory.cpp
  src/drt_structure_new/str_predict_generic.cpp
  src/drt_structure_new/str_predict_constdisvelaccpress.cpp
  src/drt_structure_new/str_predict_tangdis.cpp
  src/drt_structure_new/str_solver_factory.cpp
  src/drt_structure_new/str_timint_factory.cpp
  src/drt_structure_new/str_timint_base.cpp
  src/drt_structure_new/str_timint_basedataio.cpp
  src/drt_structure_new/str_timint_basedataio_runtime_vtk_output.cpp
  src/drt_structure_new/str_discretization_runtime_vtu_output_params.cpp
  src/drt_structure_new/str_timint_basedataio_runtime_vtp_output.cpp
  src/drt_structure_new/str_timint_basedataglobalstate.cpp
  src/drt_structure_new/str_timint_basedatasdyn.cpp
  src/drt_structure_new/str_timint_explicit.cpp
  src/drt_structure_new/str_timint_implicit.cpp
  src/drt_structure_new/str_timint_implicitbase.cpp
  src/drt_structure_new/str_timint_noxinterface.cpp
  src/drt_structure_new/str_timint_loca_continuation.cpp
  src/drt_structure_new/str_timint_locainterface.cpp
  src/drt_structure_new/str_nln_linearsystem_scaling.cpp
#  src/drt_structure_new/loca_str_constraint_arclength.cpp
  src/drt_structure_xstructure/str_timint_base_xfem.cpp
  src/drt_structure_xstructure/str_timint_implicitbase_xfem.cpp
)
set(OBJS_DRT_STRUCTURE_XSTRUCTURE
  src/drt_structure_xstructure/xstr_xstructure_state.cpp
  src/drt_structure_xstructure/xstr_xstructure_structure_state.cpp
  src/drt_structure_xstructure/xstr_multi_discretization_wrapper.cpp
  src/drt_structure_xstructure/xstr_multi_io.cpp
  src/drt_contact_xcontact/xcontact_cutwizard.cpp
  src/drt_contact_xcontact/xcontact_cutwizard_debug.cpp
  src/drt_contact_xcontact/xcontact_state_creator.cpp
)
set(OBJS_POST_DRT_S8CONVERT
  src/post_s8convert/post_drt_s8convert.cpp
)
set(OBJS_DRT_PARTICLE
  src/drt_particle/particle_timint.cpp
  src/drt_particle/particle_timint_expl.cpp
  src/drt_particle/particle_timint_centrdiff.cpp
  src/drt_particle/particle_timint_expleuler.cpp
  src/drt_particle/particle_timint_rk.cpp
  src/drt_particle/particle_timint_strategy.cpp
  src/drt_particle/particle_handler.cpp
  src/drt_particle/particle_algorithm.cpp
  src/drt_particle/particle_contact.cpp
  src/drt_particle/scatra_particle_coupling.cpp
  src/drt_particle/particle_resulttest.cpp
)
set(OBJS_DRT_PARTICLE_ALGORITHM
  src/drt_particle_algorithm/particle_sim.cpp
  src/drt_particle_algorithm/particle_algorithm.cpp
  src/drt_particle_algorithm/particle_timint.cpp
  src/drt_particle_algorithm/particle_input_generator.cpp
  src/drt_particle_algorithm/particle_gravity.cpp
  src/drt_particle_algorithm/particle_viscous_damping.cpp
  src/drt_particle_algorithm/particle_dirichlet_bc.cpp
  src/drt_particle_algorithm/particle_temperature_bc.cpp
  src/drt_particle_algorithm/particle_initial_field.cpp
  src/drt_particle_algorithm/particle_result_test.cpp
  src/drt_particle_algorithm/particle_algorithm_utils.cpp
)
set(OBJS_DRT_PARTICLE_INTERACTION
  src/drt_particle_interaction/particle_interaction_base.cpp
  src/drt_particle_interaction/particle_interaction_dem.cpp
  src/drt_particle_interaction/particle_interaction_dem_neighbor_pairs.cpp
  src/drt_particle_interaction/particle_interaction_dem_history_pairs.cpp
  src/drt_particle_interaction/particle_interaction_dem_contact.cpp
  src/drt_particle_interaction/particle_interaction_dem_contact_normal.cpp
  src/drt_particle_interaction/particle_interaction_dem_contact_tangential.cpp
  src/drt_particle_interaction/particle_interaction_dem_contact_rolling.cpp
  src/drt_particle_interaction/particle_interaction_dem_adhesion.cpp
  src/drt_particle_interaction/particle_interaction_dem_adhesion_law.cpp
  src/drt_particle_interaction/particle_interaction_dem_adhesion_surface_energy.cpp
  src/drt_particle_interaction/particle_interaction_sph.cpp
  src/drt_particle_interaction/particle_interaction_sph_kernel.cpp
  src/drt_particle_interaction/particle_interaction_sph_equationofstate.cpp
  src/drt_particle_interaction/particle_interaction_sph_equationofstate_bundle.cpp
  src/drt_particle_interaction/particle_interaction_sph_neighbor_pairs.cpp
  src/drt_particle_interaction/particle_interaction_sph_density.cpp
  src/drt_particle_interaction/particle_interaction_sph_density_correction.cpp
  src/drt_particle_interaction/particle_interaction_sph_pressure.cpp
  src/drt_particle_interaction/particle_interaction_sph_temperature.cpp
  src/drt_particle_interaction/particle_interaction_sph_momentum.cpp
  src/drt_particle_interaction/particle_interaction_sph_momentum_formulation.cpp
  src/drt_particle_interaction/particle_interaction_sph_artificialviscosity.cpp
  src/drt_particle_interaction/particle_interaction_sph_surface_tension.cpp
  src/drt_particle_interaction/particle_interaction_sph_boundary_particle.cpp
  src/drt_particle_interaction/particle_interaction_sph_heatsource.cpp
  src/drt_particle_interaction/particle_interaction_sph_phase_change.cpp
  src/drt_particle_interaction/particle_interaction_material_handler.cpp
)
set(OBJS_DRT_PARTICLE_ENGINE
  src/drt_particle_engine/particle_container.cpp
  src/drt_particle_engine/particle_container_bundle.cpp
  src/drt_particle_engine/particle_engine.cpp
  src/drt_particle_engine/particle_communication_utils.cpp
  src/drt_particle_engine/particle_runtime_vtp_writer.cpp
)
set(OBJS_DRT_PARTICLE_ENGINE_LIB
  src/drt_particle_engine/particle_enums.cpp
  src/drt_particle_engine/particle_object.cpp
)
set(OBJS_DRT_PARTICLE_WALL
  src/drt_particle_wall/particle_wall.cpp
  src/drt_particle_wall/particle_wall_datastate.cpp
  src/drt_particle_wall/particle_wall_result_test.cpp
)
set(OBJS_DRT_ACOU
  src/drt_acou/acou_dyn.cpp
  src/drt_acou/acou_impl_euler.cpp
  src/drt_acou/acou_impl.cpp
  src/drt_acou/acou_timeint.cpp
  src/drt_acou/acou_inv_resulttest.cpp
  src/drt_acou/acou_resulttest.cpp
  src/drt_acou/acou_expl.cpp
  src/drt_acou/acou_expl_worker.cpp
  src/drt_acou/acou_expl_worker_ader.cpp
  src/drt_acou/acou_expl_worker_ader_lts.cpp
  src/drt_acou/acou_expl_worker_sol.cpp
  src/drt_acou/acou_pml.cpp
  src/drt_acou/deal_tria_wrapper.cpp
  src/drt_acou/pat_imagereconstruction.cpp
  src/drt_acou/pat_utils.cpp
)
set(OBJS_DRT_ACOU_LIB
  src/drt_acou/acou_ele_boundary_calc.cpp
  src/drt_acou/acou_ele_calc.cpp
  src/drt_acou/acou_ele.cpp
  src/drt_acou/acou_ele_evaluate.cpp
  src/drt_acou/acou_ele_factory.cpp
  src/drt_acou/acou_ele_intfaces_calc.cpp
  src/drt_acou/acou_sol_ele_calc.cpp
  src/drt_acou/acou_sol_ele.cpp
)
set(OBJS_DRT_ELEMAG
  src/drt_elemag/elemag_dyn.cpp
  src/drt_elemag/elemag_timeint.cpp
  src/drt_elemag/elemag_resulttest.cpp
)
set(OBJS_DRT_ELEMAG_LIB
  src/drt_elemag/elemag_ele_boundary_calc.cpp
  src/drt_elemag/elemag_ele_calc.cpp
  src/drt_elemag/elemag_ele_evaluate.cpp
  src/drt_elemag/elemag_ele_factory.cpp
  src/drt_elemag/elemag_ele_intfaces_calc.cpp
  src/drt_elemag/elemag_ele.cpp
)
set(OBJS_LINALG
  src/linalg/linalg_ana.cpp
  src/linalg/linalg_blocksparsematrix.cpp
  src/linalg/linalg_condest.cpp
  src/linalg/linalg_downwindmatrix.cpp
  src/linalg/linalg_fixedsparsematrix.cpp
  src/linalg/linalg_krylov_projector.cpp
  src/linalg/linalg_mapextractor.cpp
  src/linalg/linalg_mlapi_operator.cpp
  src/linalg/linalg_precond.cpp
  src/linalg/linalg_projected_operator.cpp
  src/linalg/linalg_projected_precond.cpp
  src/linalg/linalg_serialdensematrix.cpp
  src/linalg/linalg_serialdensevector.cpp
  src/linalg/linalg_sparsematrix.cpp
  src/linalg/linalg_sparsematrixbase.cpp
  src/linalg/linalg_utils.cpp
  src/linalg/linalg_gauss.cpp
)
set(OBJS_SOLVER
  src/linalg/linalg_multiply.cpp
  src/linalg/linalg_solver.cpp
  src/solver/bgs2x2_operator.cpp
  src/solver/solver_cheapsimplepreconditioner.cpp
  src/solver/solver_preconditionertype.cpp
  src/solver/solver_directsolver.cpp
  src/solver/solver_krylovsolver.cpp
  src/solver/solver_pointpreconditioner.cpp
  src/solver/solver_blockpreconditioners.cpp
  src/solver/solver_krylovprojectionpreconditioner.cpp
  src/solver/solver_ifpackpreconditioner.cpp
  src/solver/solver_mlpreconditioner.cpp
  src/solver/solver_muelupreconditioner.cpp
  src/solver/solver_muelucontactpreconditioner.cpp
  src/solver/solver_muelucontactpreconditioner2.cpp
  src/solver/solver_muelucontactsppreconditioner.cpp
  src/solver/solver_muelucontactpenaltypreconditioner.cpp
  src/solver/solver_tekopreconditioner.cpp
  src/solver/solver_belossolver.cpp
  src/solver/solver_aztecsolver.cpp
  src/solver/solver_aztecsolver_projectedresidual.cpp
  src/solver/solver_stratimikossolver.cpp
  src/solver/teko/teko_baciepetraoperatorwrapper.cpp
  src/solver/muelu/MueLu_ContactAFilterFactory.cpp
  src/solver/muelu/MueLu_ContactTransferFactory.cpp
  src/solver/muelu/MueLu_ContactInfoFactory.cpp
  src/solver/muelu/MueLu_ContactASlaveDofFilterFactory.cpp
  src/solver/muelu/MueLu_ContactSPAggregationFactory.cpp
  src/solver/muelu/MueLu_MyTrilinosSmoother.cpp
  src/solver/muelu/MueLu_IterationAFactory.cpp
  src/solver/muelu/MueLu_SelectiveSaPFactory.cpp
  src/solver/muelu/MueLu_MeshtyingSPAmalgamationFactory.cpp
  src/solver/muelu/MueLu_ContactSPRepartitionInterface.cpp
  src/solver/muelu/MueLu_BaciFactoryFactory.cpp
  src/solver/solver_amgnxn_preconditioner.cpp
  src/solver/solver_amgnxn_vcycle.cpp
  src/solver/solver_amgnxn_hierarchies.cpp
  src/solver/solver_amgnxn_smoothers.cpp
  src/solver/solver_amgnxn_objects.cpp
  src/solver/solver_aztecoo_conditionnumber.cpp
)

set(OBJS_SOLVER_NONLIN_NOX
  src/solver_nonlin_nox/nox_nln_aux.cpp
  src/solver_nonlin_nox/nox_nln_floating_point_exception.cpp
  src/solver_nonlin_nox/nox_nln_globaldata.cpp
  src/solver_nonlin_nox/nox_nln_problem.cpp
  src/solver_nonlin_nox/nox_nln_linearsystem.cpp
  src/solver_nonlin_nox/nox_nln_linearsystem_factory.cpp
  src/solver_nonlin_nox/nox_nln_linearsystem_prepostoperator.cpp
  src/solver_nonlin_nox/nox_nln_group.cpp
  src/solver_nonlin_nox/nox_nln_group_prepostoperator.cpp
  src/solver_nonlin_nox/nox_nln_linesearch_factory.cpp
  src/solver_nonlin_nox/nox_nln_linesearch_prepostoperator.cpp
  src/solver_nonlin_nox/nox_nln_direction_factory.cpp
  src/solver_nonlin_nox/nox_nln_meritfunction_factory.cpp
  src/solver_nonlin_nox/nox_nln_solver_factory.cpp
  src/solver_nonlin_nox/nox_nln_solver_prepostop_generic.cpp
  src/solver_nonlin_nox/nox_nln_solver_linesearchbased.cpp
  src/solver_nonlin_nox/nox_nln_solver_ptc.cpp
  src/solver_nonlin_nox/nox_nln_constraint_group.cpp
  src/solver_nonlin_nox/nox_nln_linesearch_backtrack.cpp
  src/solver_nonlin_nox/nox_nln_meritfunction_lagrangian.cpp
  src/solver_nonlin_nox/nox_nln_meritfunction_infeasibility.cpp
  src/solver_nonlin_nox/nox_nln_direction_newton.cpp
  src/solver_nonlin_nox/nox_nln_direction_modified_newton.cpp
  src/solver_nonlin_nox/nox_nln_direction_defaultsteptest.cpp
  src/solver_nonlin_nox/nox_nln_statustest_factory.cpp
  src/solver_nonlin_nox/nox_nln_statustest_combo.cpp
  src/solver_nonlin_nox/nox_nln_statustest_normf.cpp
  src/solver_nonlin_nox/nox_nln_statustest_normupdate.cpp
  src/solver_nonlin_nox/nox_nln_statustest_normwrms.cpp
  src/solver_nonlin_nox/nox_nln_statustest_activeset.cpp
  src/solver_nonlin_nox/nox_nln_inner_statustest_factory.cpp
  src/solver_nonlin_nox/nox_nln_inner_statustest_generic.cpp
  src/solver_nonlin_nox/nox_nln_inner_statustest_armijo.cpp
  src/solver_nonlin_nox/nox_nln_inner_statustest_filter.cpp
  src/solver_nonlin_nox/nox_nln_inner_statustest_upperbound.cpp
  src/solver_nonlin_nox/nox_nln_inner_statustest_combo.cpp
  src/solver_nonlin_nox/nox_nln_inner_statustest_volume_change.cpp
  src/drt_structure_new/nox_nln_str_linearsystem.cpp
  src/drt_contact_aug/nox_nln_contact_linearsystem.cpp
  src/drt_contact/nox_nln_meshtying_linearsystem.cpp
  src/drt_cardiovascular0d/nox_nln_cardiovascular0d_linearsystem.cpp
  src/drt_constraint/nox_nln_lagpenconstraint_linearsystem.cpp
  src/drt_scatra/nox_nln_scatra_linearsystem.cpp
)

set(OBJS_LOCA_CONTINUATION
  src/loca_continuation/loca_nln_problem.cpp
  src/loca_continuation/loca_nln_group.cpp
  src/loca_continuation/loca_nln_stepper.cpp
  src/loca_continuation/loca_nln_statustest_factory.cpp
)

set(OBJS_DRT_CARDIOVASCULAR0D
  src/drt_cardiovascular0d/cardiovascular0d.cpp
  src/drt_cardiovascular0d/cardiovascular0d_manager.cpp
  src/drt_cardiovascular0d/cardiovascular0d_dofset.cpp
  src/drt_cardiovascular0d/cardiovascular0d_resulttest.cpp
  src/drt_cardiovascular0d/cardiovascular0d_4elementwindkessel.cpp
  src/drt_cardiovascular0d/cardiovascular0d_arterialproxdist.cpp
  src/drt_cardiovascular0d/cardiovascular0d_syspulcirculation.cpp
  src/drt_cardiovascular0d/cardiovascularrespiratory0d_syspulperiphcirculation.cpp
)

set(OBJS_DRT_CONSTRAINT
  src/drt_constraint/constraint.cpp
  src/drt_constraint/constraint_manager.cpp
  src/drt_constraint/constraintdofset.cpp
  src/drt_constraint/constraintpenalty.cpp
  src/drt_constraint/constraintsolver.cpp
  src/drt_constraint/monitor.cpp
  src/drt_constraint/multipointconstraint.cpp
  src/drt_constraint/multipointconstraint2.cpp
  src/drt_constraint/multipointconstraint3.cpp
  src/drt_constraint/multipointconstraint3penalty.cpp
  src/drt_constraint/lagpenconstraint_noxinterface.cpp
  src/drt_constraint/springdashpot.cpp
  src/drt_constraint/springdashpot_manager.cpp
  src/drt_constraint/springdashpot_new.cpp
)

set(OBJS_DRT_IO_LIB
  src/drt_io/io.cpp
  src/drt_io/io_control.cpp
  src/drt_io/io_hdf.cpp
  src/drt_io/io_gmsh.cpp
  src/drt_io/every_iteration_writer.cpp
  src/drt_io/vtk_writer_base.cpp
  src/drt_io/vtu_writer.cpp
  src/drt_io/vtp_writer.cpp
  src/drt_io/runtime_vtu_writer.cpp
  src/drt_io/runtime_vtp_writer.cpp
  src/drt_io/discretization_runtime_vtu_writer.cpp
  src/drt_io/discretization_runtime_vtp_writer.cpp
  src/drt_beam3/beam_discretization_runtime_vtu_writer.cpp
  src/drt_io/io_pstream.cpp
  src/drt_lib/drt_domainreader.cpp
  src/drt_lib/drt_elementreader.cpp
  src/drt_lib/drt_inputreader.cpp
  src/drt_lib/drt_linedefinition.cpp
  src/drt_lib/drt_materialdefinition.cpp
  src/drt_lib/drt_nodereader.cpp
  src/drt_lib/drt_particlereader.cpp
)
set(OBJS_FORTRAN
  src/fortran/hex_fserv.f
)
set(OBJS_DRT_NURBS_DISCRET
  src/drt_nurbs_discret/drt_apply_nurbs_initial_condition.cpp
)
set(OBJS_DRT_B3_LIB
  src/drt_beam3/beam_discretization_runtime_vtu_output_params.cpp
  src/drt_beam3/triad_interpolation.cpp
  src/drt_beam3/triad_interpolation_local_rotation_vectors.cpp
  src/drt_beam3/beam3_base.cpp
  src/drt_beam3/beam3.cpp
  src/drt_beam3/beam3_evaluate.cpp
  src/drt_beam3/beam3_input.cpp
  src/drt_beam3/beam3r.cpp
  src/drt_beam3/beam3r_evaluate.cpp
  src/drt_beam3/beam3r_input.cpp
  src/drt_beam3/beam3eb.cpp
  src/drt_beam3/beam3eb_evaluate.cpp
  src/drt_beam3/beam3eb_input.cpp
  src/drt_beam3/beam3k.cpp
  src/drt_beam3/beam3k_evaluate.cpp
  src/drt_beam3/beam3k_input.cpp
)
set(OBJS_DRT_FLUID_ELE_LIB
  src/drt_fluid_ele/fluid_ele.cpp
  src/drt_fluid_ele/fluid_ele_immersed.cpp
  src/drt_fluid_ele/fluid_ele_poro_immersed.cpp
  src/drt_fluid_ele/fluid_ele_immersed_base.cpp
  src/drt_fluid_ele/fluid_ele_boundary_calc.cpp
  src/drt_fluid_ele/fluid_ele_boundary_calc_std.cpp
  src/drt_fluid_ele/fluid_ele_boundary_calc_poro.cpp
  src/drt_fluid_ele/fluid_ele_boundary_parent_calc.cpp
  src/drt_fluid_ele/fluid_ele_boundary.cpp
  src/drt_fluid_ele/fluid_ele_boundary_evaluate.cpp
  src/drt_fluid_ele/fluid_ele_poro_boundary.cpp
  src/drt_fluid_ele/fluid_ele_poro_boundary_evaluate.cpp
  src/drt_fluid_ele/fluid_ele_calc.cpp
  src/drt_fluid_ele/fluid_ele_calc_newOST.cpp
  src/drt_fluid_ele/fluid_ele_calc_xwall_quadrature_service.cpp
  src/drt_fluid_ele/fluid_ele_calc_xwall.cpp
  src/drt_fluid_ele/fluid_ele_xwall.cpp
  src/drt_fluid_ele/fluid_ele_xwall_boundary.cpp
  src/drt_fluid_ele/fluid_ele_calc_ale_service.cpp
  src/drt_fluid_ele/fluid_ele_calc_general_service.cpp
  src/drt_fluid_ele/fluid_ele_calc_loma.cpp
  src/drt_fluid_ele/fluid_ele_intfaces.cpp
  src/drt_fluid_ele/fluid_ele_calc_intfaces_stab.cpp
  src/drt_fluid_ele/fluid_ele_intfaces_calc.cpp
  src/drt_fluid_ele/fluid_ele_calc_loma_service.cpp
  src/drt_fluid_ele/fluid_ele_calc_topopt_service.cpp
  src/drt_fluid_ele/fluid_ele_calc_poro.cpp
  src/drt_fluid_ele/fluid_ele_calc_poro_p1.cpp
  src/drt_fluid_ele/fluid_ele_calc_poro_p1_immersed.cpp
  src/drt_fluid_ele/fluid_ele_calc_immersed.cpp
  src/drt_fluid_ele/fluid_ele_calc_std.cpp
  src/drt_fluid_ele/fluid_ele_calc_turbulence_service.cpp
  src/drt_fluid_ele/fluid_ele_calc_xfem.cpp
  src/drt_fluid_ele/fluid_ele_calc_xfem_coupling.cpp
  src/drt_fluid_ele/fluid_ele_calc_xfem_coupling_impl.cpp
  src/drt_fluid_ele/fluid_ele_evaluate.cpp
  src/drt_fluid_ele/fluid_ele_factory.cpp
  src/drt_fluid_ele/fluid_ele_boundary_factory.cpp
  src/drt_fluid_ele/fluid_ele_input.cpp
  src/drt_fluid_ele/fluid_ele_parameter.cpp
  src/drt_fluid_ele/fluid_ele_parameter_std.cpp
  src/drt_fluid_ele/fluid_ele_parameter_poro.cpp
  src/drt_fluid_ele/fluid_ele_parameter_xfem.cpp
  src/drt_fluid_ele/fluid_ele_parameter_timint.cpp
  src/drt_fluid_ele/fluid_ele_parameter_intface.cpp
  src/drt_fluid_ele/fluid_ele_tds.cpp
  src/drt_fluid_ele/fluid_ele_poro.cpp
  src/drt_fluid_ele/fluid_ele_poro_evaluate.cpp
  src/drt_fluid_ele/fluid_ele_hdg.cpp
  src/drt_fluid_ele/fluid_ele_calc_hdg.cpp
  src/drt_opti/topopt_fluidAdjoint3_boundary.cpp
  src/drt_opti/topopt_fluidAdjoint3_impl.cpp
  src/drt_opti/topopt_fluidAdjoint3_impl_parameter.cpp
)
set(OBJS_DRT_T3_LIB
  src/drt_truss3/truss3.cpp
  src/drt_truss3/truss3_evaluate.cpp
  src/drt_truss3/truss3_input.cpp
)
set(OBJS_DRT_DISCSH3_LIB
  src/drt_discsh3/discsh3.cpp
  src/drt_discsh3/discsh3_line.cpp
  src/drt_discsh3/discsh3_evaluate.cpp
  src/drt_discsh3/discsh3_line_evaluate.cpp
  src/drt_discsh3/discsh3_input.cpp
)
set(OBJS_DRT_SP3_LIB
  src/drt_spring3/spring3.cpp
  src/drt_spring3/spring3_evaluate.cpp
  src/drt_spring3/spring3_input.cpp
)
set(OBJS_DRT_T3CL_LIB
  src/drt_truss3cl/truss3cl.cpp
  src/drt_truss3cl/truss3cl_evaluate.cpp
  src/drt_truss3cl/truss3cl_input.cpp
)
set(OBJS_DRT_ALEALG
  src/drt_ale/ale.cpp
  src/drt_ale/ale_dyn.cpp
  src/drt_ale/ale_meshsliding.cpp
  src/drt_ale/ale_meshtying.cpp
  src/drt_ale/ale_resulttest.cpp
  src/drt_ale/ale_utils_mapextractor.cpp
  src/drt_ale/ale_utils_clonestrategy.cpp
)
set(OBJS_DRT_W1_LIB
  src/drt_w1/wall1.cpp
  src/drt_w1/wall1_eas.cpp
  src/drt_w1/wall1_evaluate.cpp
  src/drt_w1/wall1_gemm.cpp
  src/drt_w1/wall1_input.cpp
  src/drt_w1/wall1_line.cpp
  src/drt_w1/wall1_line_evaluate.cpp
  src/drt_w1/wall1_mat.cpp
  src/drt_w1/wall1_nurbs.cpp
  src/drt_w1/wall1_poro.cpp
  src/drt_w1/wall1_poro_eletypes.cpp
  src/drt_w1/wall1_poro_evaluate.cpp
  src/drt_w1/wall1_poro_p1.cpp
  src/drt_w1/wall1_poro_p1_eletypes.cpp
  src/drt_w1/wall1_poro_p1_evaluate.cpp
  src/drt_w1/wall1_poro_p1_scatra.cpp
  src/drt_w1/wall1_poro_p1_scatra_eletypes.cpp
  src/drt_w1/wall1_poro_scatra.cpp
  src/drt_w1/wall1_poro_scatra_eletypes.cpp
  src/drt_w1/wall1_scatra.cpp
  src/drt_w1/wall1_scatra_evaluate.cpp
)
set(OBJS_DRT_CONTACT
  src/drt_contact/contact_abstract_strategy.cpp
  src/drt_contact/contact_noxinterface.cpp
  src/drt_contact/contact_coupling2d.cpp
  src/drt_contact/contact_coupling3d.cpp
  src/drt_contact/contact_line_coupling.cpp
  src/drt_contact/contact_integrator.cpp
  src/drt_contact/contact_nitsche_integrator.cpp
  src/drt_contact/contact_nitsche_integrator_tsi.cpp
  src/drt_contact/contact_nitsche_integrator_poro.cpp
  src/drt_contact/contact_ehl_integrator.cpp
  src/drt_contact/contact_integrator_factory.cpp
  src/drt_contact/contact_integrator_assemble.cpp
  src/drt_contact/contact_interface.cpp
  src/drt_contact/contact_interface_tools.cpp
  src/drt_contact/contact_interpolator.cpp
  src/drt_contact/contact_lagrange_strategy.cpp
  src/drt_contact/contact_monocoupled_lagrange_strategy.cpp
  src/drt_contact/contact_poro_lagrange_strategy.cpp
  src/drt_contact/contact_tsi_lagrange_strategy.cpp
  src/drt_contact/contact_tsi_interface.cpp
  src/drt_contact/contact_manager.cpp
  src/drt_contact/contact_strategy_factory.cpp
  src/drt_contact/contact_penalty_strategy.cpp
  src/drt_contact/contact_nitsche_strategy.cpp
  src/drt_contact/contact_nitsche_strategy_tsi.cpp
  src/drt_contact/contact_nitsche_strategy_poro.cpp
  src/drt_contact/contact_utils.cpp
  src/drt_contact/contact_wear_interface.cpp
  src/drt_contact/contact_wear_interface_tools.cpp
  src/drt_contact/contact_wear_lagrange_strategy.cpp
  src/drt_contact/meshtying_strategy_factory.cpp
  src/drt_contact/meshtying_abstract_strategy.cpp
  src/drt_contact/meshtying_noxinterface.cpp
  src/drt_contact/meshtying_contact_bridge.cpp
  src/drt_contact/meshtying_lagrange_strategy.cpp
  src/drt_contact/meshtying_poro_lagrange_strategy.cpp
  src/drt_contact/meshtying_manager.cpp
  src/drt_contact/meshtying_penalty_strategy.cpp
  src/drt_contact/selfcontact_binarytree.cpp
  src/drt_contact/selfcontact_binarytree_unbiased.cpp
)
set(OBJS_DRT_CONTACT_AUG
  src/drt_contact_aug/contact_augmented_strategy.cpp
  src/drt_contact_aug/contact_augmented_strategy_tools.cpp
  src/drt_contact_aug/contact_augmented_interface.cpp
  src/drt_contact_aug/contact_augmented_interface_tools.cpp
  src/drt_contact_aug/contact_augmented_integrator.cpp
  src/drt_contact_aug/contact_augmented_integrator_evaluate.cpp
  src/drt_contact_aug/contact_aug_active_set.cpp
  src/drt_contact_aug/contact_aug_projector.cpp
  src/drt_contact_aug/contact_aug_lagrange_strategy.cpp
  src/drt_contact_aug/contact_aug_lagrange_interface.cpp
  src/drt_contact_aug/contact_aug_steepest_ascent_sp_strategy.cpp
  src/drt_contact_aug/contact_aug_steepest_ascent_strategy.cpp
  src/drt_contact_aug/contact_aug_steepest_ascent_interface.cpp
  src/drt_contact_aug/contact_aug_penalty_update.cpp
  src/drt_contact_aug/contact_aug_combo_strategy.cpp
  src/drt_contact_aug/contact_aug_combo_strategy_switching.cpp
  src/drt_contact_aug/contact_aug_potential.cpp
  src/drt_contact_aug/contact_aug_plot.cpp
  src/drt_contact_aug/contact_integrator_utils.cpp
  src/drt_contact_aug/contact_aug_integrator_policy.cpp
  src/drt_contact_aug/contact_aug_integrator_policy_debug.cpp
  src/drt_contact_aug/contact_aug_lagrange_multiplier_function.cpp
  src/drt_contact_aug/contact_aug_interface_assembly.cpp
  src/drt_contact_aug/contact_aug_timemonitor.cpp
  src/drt_contact_aug/contact_aug_parallel_distribution_controller.cpp
  src/drt_contact_xcontact/xcontact_strategy.cpp
  src/drt_contact_xcontact/xcontact_interface.cpp
  src/drt_contact_xcontact/xcontact_integrator.cpp
)
set(OBJS_DRT_CUT
  src/drt_cut/cut_options.cpp
  src/drt_cut/base_vol.cpp
  src/drt_cut/cut_boundarycell.cpp
  src/drt_cut/cut_boundingbox.cpp
  src/drt_cut/cut_coloredgraph.cpp
  src/drt_cut/cut_cycle.cpp
  src/drt_cut/direct_divergence.cpp
  src/drt_cut/direct_divergence_refplane.cpp
  src/drt_cut/cut_edge.cpp
  src/drt_cut/cut_element.cpp
  src/drt_cut/cut_position.cpp
  src/drt_cut/cut_intersection.cpp
  src/drt_cut/cut_elementhandle.cpp
  src/drt_cut/cut_facet.cpp
  src/drt_cut/cut_facetgraph.cpp
  src/drt_cut/cut_find_cycles.cpp
  src/drt_cut/cut_integrationcell.cpp
  src/drt_cut/cut_integrationcellcreator.cpp
  src/drt_cut/cut_kernel.cpp
  src/drt_cut/cut_combintersection.cpp
  src/drt_cut/cut_levelsetintersection.cpp
  src/drt_cut/cut_levelsetside.cpp
  src/drt_cut/cut_line.cpp
  src/drt_cut/cut_mesh.cpp
  src/drt_cut/cut_meshhandle.cpp
  src/drt_cut/cut_meshintersection.cpp
  src/drt_cut/cut_node.cpp
  src/drt_cut/cut_output.cpp
  src/drt_cut/cut_parallel.cpp
  src/drt_cut/cut_parentintersection.cpp
  src/drt_cut/cut_point.cpp
  src/drt_cut/cut_pointgraph.cpp
  src/drt_cut/cut_pointpool.cpp
  src/drt_cut/cut_selfcut.cpp
  src/drt_cut/cut_side.cpp
  src/drt_cut/cut_sidehandle.cpp
  src/drt_cut/cut_tetmesh.cpp
  src/drt_cut/cut_tetmeshintersection.cpp
  src/drt_cut/cut_triangulateFacet.cpp
  src/drt_cut/cut_volumecell.cpp
  src/drt_cut/line_integration.cpp
  src/drt_cut/facet_integration.cpp
  src/drt_cut/boundarycell_integration.cpp
  src/drt_cut/volume_integration.cpp
  src/drt_cut/least_squares.cpp
  src/drt_cut/quadrature_compression.cpp
  src/drt_cut/cut_cutwizard.cpp
  src/drt_cut/cut_pointgraph_simple.cpp
  src/drt_cut/cut_facetgraph_simple.cpp
  src/drt_cut/cut_clnwrapper.cpp
  src/drt_cut/cut_memory_manager.cpp
)
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
  tests/cut_test/drt_singletondestruction.cpp
  tests/cut_test/cut_test_intersection.cpp
  tests/cut_test/drt_singletondestruction.cpp
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
  src/drt_fem_general/drt_utils_gausspoints.cpp
  src/drt_fem_general/drt_utils_fem_shapefunctions.cpp
  src/drt_fem_general/drt_utils_integration.cpp
  src/drt_cut/cut_clnwrapper.cpp
  src/drt_cut/cut_memory_manager.cpp 
)
set(OBJS_GRAPH_TEST
  tests/graph_test/graph_test_main.cpp
  tests/graph_test/graph_test_read.cpp
  tests/graph_test/graph_test_gnuplot.cpp
  tests/graph_test/graph_test_output.cpp
  tests/graph_test/graph_test_find_cycles.cpp
)
set(OBJS_DRT_WEAR
  src/drt_wear/wear_algorithm.cpp
  src/drt_wear/wear_dyn.cpp
  src/drt_wear/wear_partitioned.cpp
)
set(OBJS_DRT_POROELAST
  src/drt_poroelast/poro_dyn.cpp
  src/drt_poroelast/poro_base.cpp
  src/drt_poroelast/poroelast_monolithic.cpp
  src/drt_poroelast/poroelast_utils.cpp
  src/drt_poroelast/poro_utils_clonestrategy.cpp
  src/drt_poroelast/poro_partitioned.cpp
  src/drt_poroelast/poro_monolithicsplit.cpp
  src/drt_poroelast/poro_monolithicstructuresplit.cpp
  src/drt_poroelast/poro_monolithicfluidsplit.cpp
  src/drt_poroelast/poro_monolithicsplit_nopenetration.cpp
  src/drt_poroelast/poro_monolithicmeshtying.cpp
  src/drt_poroelast/poro_scatra_base.cpp
  src/drt_poroelast/poro_scatra_part.cpp
  src/drt_poroelast/poro_scatra_part_1wc.cpp
  src/drt_poroelast/poro_scatra_part_2wc.cpp
  src/drt_poroelast/poro_scatra_monolithic.cpp
)

set(OBJS_DRT_POROFLUIDMULTIPHASE
  src/drt_porofluidmultiphase/porofluidmultiphase_dyn.cpp
  src/drt_porofluidmultiphase/porofluidmultiphase_timint_implicit.cpp
  src/drt_porofluidmultiphase/porofluidmultiphase_timint_ost.cpp
  src/drt_porofluidmultiphase/porofluidmultiphase_meshtying_strategy_std.cpp
  src/drt_porofluidmultiphase/porofluidmultiphase_meshtying_strategy_artery.cpp
  src/drt_porofluidmultiphase/porofluidmultiphase_resulttest.cpp
  src/drt_porofluidmultiphase/porofluidmultiphase_utils.cpp
)

set(OBJS_DRT_POROMULTIPHASE
  src/drt_poromultiphase/poromultiphase_dyn.cpp
  src/drt_poromultiphase/poromultiphase_base.cpp
  src/drt_poromultiphase/poromultiphase_partitioned.cpp
  src/drt_poromultiphase/poromultiphase_monolithic.cpp
  src/drt_poromultiphase/poromultiphase_partitioned_twoway.cpp
  src/drt_poromultiphase/poromultiphase_monolithic_twoway.cpp
  src/drt_poromultiphase/poromultiphase_utils_clonestrategy.cpp
  src/drt_poromultiphase/poromultiphase_utils.cpp
  src/drt_poromultiphase_scatra/poromultiphase_scatra_dyn.cpp
  src/drt_poromultiphase_scatra/poromultiphase_scatra_utils.cpp
  src/drt_poromultiphase_scatra/poromultiphase_scatra_base.cpp
  src/drt_poromultiphase_scatra/poromultiphase_scatra_partitioned_twoway.cpp
  src/drt_poromultiphase_scatra/poromultiphase_scatra_monolithic_twoway.cpp
)

set(OBJS_DRT_POROMULTIPHASE_ART_COUPL
  src/drt_poromultiphase_scatra/poromultiphase_scatra_artery_coupling_base.cpp
  src/drt_poromultiphase_scatra/poromultiphase_scatra_artery_coupling_nodebased.cpp
  src/drt_poromultiphase_scatra/poromultiphase_scatra_artery_coupling_linebased.cpp
  src/drt_poromultiphase_scatra/poromultiphase_scatra_artery_coupling_pair.cpp
)

set(OBJS_DRT_FPSI
  src/drt_fpsi/fpsi_utils.cpp
  src/drt_fpsi/fpsi_dyn.cpp
  src/drt_fpsi/fpsi_monolithic.cpp
  src/drt_fpsi/fpsi_monolithic_plain.cpp
  src/drt_fpsi/fpsi.cpp
  src/drt_fpsi/fpsi_coupling.cpp
)

set(OBJS_DRT_SSI
src/drt_ssi/ssi_dyn.cpp
src/drt_ssi/ssi_base.cpp
src/drt_ssi/ssi_clonestrategy.cpp
src/drt_ssi/ssi_coupling.cpp
src/drt_ssi/ssi_monolithic.cpp
src/drt_ssi/ssi_monolithic_convcheck_strategies.cpp
src/drt_ssi/ssi_partitioned.cpp
src/drt_ssi/ssi_partitioned_1wc.cpp
src/drt_ssi/ssi_partitioned_2wc.cpp
src/drt_ssi/ssi_resulttest.cpp
src/drt_ssi/ssi_utils.cpp
src/drt_ssi/ssi_str_model_evaluator_monolithic.cpp
src/drt_ssi/ssi_str_model_evaluator_partitioned.cpp
)

set(OBJS_DRT_IMMERSED_LIB
src/drt_fsi/fsi_nox_aitken_immersed_ale.cpp
src/drt_immersed_problem/immersed_node.cpp
src/drt_immersed_problem/immersed_discretization.cpp
src/drt_immersed_problem/str_timint_ost_immersed.cpp
src/drt_immersed_problem/immersed_field_exchange_manager.cpp
src/drt_immersed_problem/scatra_timint_ost_endoexocytosis.cpp
)

set(OBJS_DRT_IMMERSED_PROBLEM
src/drt_immersed_problem/immersed_base.cpp
src/drt_immersed_problem/immersed_problem_dyn.cpp
src/drt_immersed_problem/immersed_partitioned.cpp
src/drt_immersed_problem/immersed_partitioned_fsi.cpp
src/drt_immersed_problem/immersed_partitioned_multiphysics.cpp
src/drt_immersed_problem/immersed_partitioned_confine_cell.cpp
src/drt_immersed_problem/immersed_partitioned_cellmigration.cpp
src/drt_immersed_problem/immersed_partitioned_cellcontraction.cpp
src/drt_immersed_problem/immersed_partitioned_adhesion_traction.cpp
src/drt_immersed_problem/immersed_partitioned_protrusion_formation.cpp
src/drt_immersed_problem/immersed_partitioned_fsi_dirichletneumann.cpp
src/drt_immersed_problem/immersed_partitioned_flow_cell_interaction.cpp
src/drt_immersed_problem/immersed_partitioned_fsi_dirichletneumann_ale.cpp
src/drt_immersed_problem/immersed_partitioned_fsi_dirichletneumann_membrane.cpp
src/drt_immersed_problem/ssi_partitioned_2wc_protrusionformation.cpp
src/drt_immersed_problem/ssi_partitioned_2wc_biochemomechano.cpp
src/drt_immersed_problem/ssi_partitioned_2wc_adhesiondynamics.cpp
src/drt_immersed_problem/fsi_partitioned_immersed.cpp
)

set(OBJS_DRT_MOR
  src/drt_mor/mor_pod.cpp
)

set(OBJS_DRT_PASI
src/drt_pasi/pasi_dyn.cpp
src/drt_pasi/pasi_partitioned.cpp
src/drt_pasi/pasi_partitioned_onewaycoup.cpp
src/drt_pasi/pasi_partitioned_twowaycoup.cpp
src/drt_pasi/pasi_utils.cpp
src/drt_pasi/pasi_str_model_evaluator_partitioned.cpp
)

set(OBJS_DRT_TUTORIAL
src/drt_tutorial/tutorial_dyn.cpp
src/drt_tutorial/tutorial_utils.cpp
src/drt_tutorial/tutorial_nln_truss.cpp
src/drt_tutorial/tutorial_struct_ale_coupling.cpp
src/drt_tutorial/tutorial_fixedpoint_scheme.cpp
src/drt_tutorial/tutorial_material/tutorial_mat_stvenantkirchhoff.cpp
src/drt_tutorial/tutorial_ele/tutorial_ele_calc.cpp
src/drt_tutorial/tutorial_ele/tutorial_ele.cpp
)

set(OBJS_DRT_FIBER_LIB
src/drt_fiber/drt_fiber_node.cpp
)

set(OBJS_DRT_GEOMETRY_PAIR
src/drt_geometry_pair/geometry_pair_factory.cpp
src/drt_geometry_pair/geometry_pair.cpp
src/drt_geometry_pair/geometry_pair_utility_classes.cpp
src/drt_geometry_pair/geometry_pair_utility_functions.cpp
src/drt_geometry_pair/geometry_pair_evaluation_data_global.cpp
src/drt_geometry_pair/geometry_pair_evaluation_data_base.cpp
src/drt_geometry_pair/geometry_pair_line_to_volume_evaluation_data.cpp
src/drt_geometry_pair/geometry_pair_line_to_volume.cpp
src/drt_geometry_pair/geometry_pair_line_to_volume_gauss_point_projection.cpp
src/drt_geometry_pair/geometry_pair_line_to_volume_gauss_point_projection_cross_section.cpp
src/drt_geometry_pair/geometry_pair_line_to_volume_segmentation.cpp
)
set(OBJS_DRT_MIXTURE
  src/drt_mixture/mixture_constituent.cpp
  src/drt_mixture/mixture_rule.cpp
  src/drt_mixture/mixture_constituent_elasthyper.cpp)


###------------------------------------------------------------------ Element flags

#   all possible elements

option (D_ALE         "ale element (2D+3D)" ON)
if (D_ALE)
  add_definitions("-DD_ALE")
endif(D_ALE)

option (D_FLUID3      "fluid3" ON)
if (D_FLUID3)
  add_definitions("-DD_FLUID3")
endif(D_FLUID3)

option (D_SHELL8      "shell8" ON)
if (D_SHELL8)
  add_definitions("-DD_SHELL8")
endif(D_SHELL8)

option (D_THERMO      "thermal elements" ON)
if (D_THERMO)
  add_definitions("-DD_THERMO")
endif(D_THERMO)

option (D_ALE_BFLOW   "ALE deformation included in flow calculation" OFF)
if (D_ALE_BFLOW)
  add_definitions("-DD_ALE_BFLOW")
endif(D_ALE_BFLOW)


###------------------------------------------------------------------ Binaries

# hide libraries

if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/src)
endif (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

# Define the libraries. The link order is important!

# Compile element libraries first. There are files that take time to compile.
add_library(ele_fluid ${OBJS_DRT_FLUID_ELE_LIB})
add_library(ele_ale ${OBJS_DRT_ALE_LIB})
add_library(ele_b3 ${OBJS_DRT_B3_LIB})
add_library(ele_rigidsphere ${OBJS_DRT_RIGIDSPHERE_LIB})
add_library(ele_bele3 ${OBJS_DRT_BELE3_LIB})
add_library(ele_opti ${OBJS_DRT_TOPOPT_LIB})
add_library(ele_s8 ${OBJS_DRT_S8_LIB})
add_library(ele_membrane ${OBJS_DRT_MEMBRANE_LIB})
add_library(ele_lubrication ${OBJS_DRT_LUBRICATION_LIB})
add_library(ele_porofluidmultiphase ${OBJS_DRT_POROFLUIDMULTIPHASE_LIB})
add_library(ele_scatra ${OBJS_DRT_SCATRA_LIB})
add_library(ele_discsh3 ${OBJS_DRT_DISCSH3_LIB})
add_library(ele_sp3 ${OBJS_DRT_SP3_LIB})
add_library(ele_so ${OBJS_SOH8_LIB})
add_library(ele_t3 ${OBJS_DRT_T3_LIB})
add_library(ele_t3cl ${OBJS_DRT_T3CL_LIB})
add_library(ele_thermo ${OBJS_DRT_THERMO_LIB})
add_library(ele_to3 ${OBJS_DRT_TO3_LIB})
add_library(ele_w1 ${OBJS_DRT_W1_LIB})
add_library(ele_acou ${OBJS_DRT_ACOU_LIB})
add_library(ele_elemag ${OBJS_DRT_ELEMAG_LIB})

set(ELEMENT_LIBRARIES
  ele_ale
  ele_b3
  ele_rigidsphere
  ele_bele3
  ele_discsh3
  ele_fluid
  ele_opti
  ele_s8
  ele_membrane
  ele_lubrication
  ele_porofluidmultiphase
  ele_scatra
  ele_sp3
  ele_so
  ele_t3
  ele_t3cl
  ele_thermo
  ele_to3
  ele_w1
  ele_acou
  ele_elemag
  )

# algorithm libraries that depend on adapters
add_library(alg_elch ${OBJS_DRT_ELCH})
add_library(alg_var ${OBJS_DRT_VAR})
add_library(alg_fsi ${OBJS_DRT_FSI})
add_library(alg_immersed_problem ${OBJS_DRT_IMMERSED_PROBLEM})
add_library(alg_fsi_xfem ${OBJS_DRT_FSI_XFEM})
add_library(alg_loma ${OBJS_DRT_LOMA})
add_library(alg_cavitation ${OBJS_DRT_CAVITATION})
add_library(alg_two_phase ${OBJS_DRT_TWO_PHASE})
add_library(alg_contact_xfem ${OBJS_DRT_CONTACT_XCONTACT})
add_library(alg_tsi ${OBJS_DRT_TSI})
add_library(alg_sti ${OBJS_DRT_STI})
add_library(alg_wear ${OBJS_DRT_WEAR})
add_library(alg_poroelast ${OBJS_DRT_POROELAST})
add_library(alg_poromultiphase ${OBJS_DRT_POROMULTIPHASE})
add_library(alg_fpsi ${OBJS_DRT_FPSI})
add_library(alg_fs3i ${OBJS_DRT_FS3I})
add_library(alg_topopt ${OBJS_DRT_TOPOPT})
add_library(alg_ssi ${OBJS_DRT_SSI})
add_library(alg_ehl ${OBJS_DRT_EHL})
add_library(alg_pasi ${OBJS_DRT_PASI})

set(MULTIFIELD_LIBRARIES
  alg_elch
  alg_var
  alg_fsi
  alg_fsi_xfem
  alg_loma
  alg_cavitation
  alg_two_phase
  alg_contact_xfem
  alg_tsi
  alg_sti
  alg_wear
  alg_poroelast
  alg_poromultiphase
  alg_immersed_problem
  alg_fpsi
  alg_fs3i
  alg_topopt
  alg_ssi
  alg_ehl
  alg_pasi
  )

# the adapter layer
add_library(drt_adapter ${OBJS_DRT_ADAPT})

# algorithm libraries the adapters depend uo
add_library(alg_ale ${OBJS_DRT_ALEALG})
add_library(alg_artnet ${OBJS_DRT_ARTNET})
#add_library(alg_contact SHARED ${OBJS_DRT_CONTACT})
add_library(alg_contact ${OBJS_DRT_CONTACT})
add_library(alg_contact_aug ${OBJS_DRT_CONTACT_AUG})
add_library(alg_fluid ${OBJS_DRT_FLUID})
add_library(alg_fluid_turbulence ${OBJS_DRT_FLUID_TURBULENCE})
add_library(alg_xfluid ${OBJS_DRT_FLUID_XFLUID})
add_library(alg_invanalysis ${OBJS_INV_ANALYSIS})
add_library(alg_statinvanalysis ${OBJS_STATINV_ANALYSIS})
add_library(alg_mlmc ${OBJS_MLMC})
add_library(alg_redairways ${OBJS_DRT_RED_AIRWAYS})
add_library(alg_lubrication ${OBJS_DRT_LUBRICATION})
add_library(alg_porofluidmultiphase ${OBJS_DRT_POROFLUIDMULTIPHASE})
add_library(alg_scatra ${OBJS_DRT_SCATRA})
add_library(alg_structure ${OBJS_DRT_STRUCTURE})
add_library(alg_xstructure ${OBJS_DRT_STRUCTURE_XSTRUCTURE})
add_library(alg_thermo ${OBJS_DRT_THERMO})
add_library(alg_fluid_adjoint ${OBJS_DRT_TOPOPT_ADJOINT})
add_library(alg_opti ${OBJS_OPTI})
add_library(alg_particle ${OBJS_DRT_PARTICLE})
add_library(alg_particle_algorithm ${OBJS_DRT_PARTICLE_ALGORITHM})
add_library(alg_particle_interaction ${OBJS_DRT_PARTICLE_INTERACTION})
add_library(alg_particle_engine ${OBJS_DRT_PARTICLE_ENGINE})
add_library(alg_particle_wall ${OBJS_DRT_PARTICLE_WALL})
add_library(alg_acou ${OBJS_DRT_ACOU})
add_library(alg_elemag ${OBJS_DRT_ELEMAG})
add_library(alg_tutorial ${OBJS_DRT_TUTORIAL})

set(SINGLEFIELD_LIBRARIES
  alg_ale
  alg_artnet
  alg_fluid
  alg_fluid_turbulence
  alg_fluid
  alg_xfluid
  alg_redairways
  alg_lubrication
  alg_porofluidmultiphase
  alg_scatra
  alg_structure
  alg_xstructure
  alg_invanalysis
  alg_statinvanalysis
  alg_mlmc
  alg_thermo
  alg_contact
  alg_contact_aug
  alg_fluid_adjoint
  alg_opti
  alg_particle
  alg_particle_algorithm
  alg_particle_interaction
  alg_particle_engine
  alg_particle_wall
  alg_acou
  alg_elemag
  alg_tutorial
  )

# general framework
add_library(drt_nurbs ${OBJS_DRT_NURBS_DISCRET})
add_library(drt_constraint ${OBJS_DRT_CONSTRAINT})
add_library(drt_cardiovascular0d ${OBJS_DRT_CARDIOVASCULAR0D})
#add_library(drt_mortar SHARED ${OBJS_DRT_MORTAR})
add_library(drt_mortar ${OBJS_DRT_MORTAR})
add_library(drt_volmortar ${OBJS_DRT_VOLMORTAR})
add_library(drt_mor ${OBJS_DRT_MOR})

# material definitions
add_library(drt_mat ${OBJS_DRT_MAT_LIB} ${OBJS_DRT_MIXTURE})
add_library(drt_micromat ${OBJS_DRT_MAT} ${OBJS_DRT_STRU_MULTI})
add_library(drt_scatra_multi ${OBJS_DRT_MAT_SCATRA_MULTI})

# special library that demands all parobject types
add_library(drt_lib_register ${OBJS_DRT_REGISTER})

# libraries that are used by baci and by filters
add_library(drt_comm ${OBJS_DRT_COMM})
add_library(drt_lib ${OBJS_DRT_LIB})
add_library(drt_lib_inpar ${OBJS_DRT_INPAR})
add_library(drt_lib_constraint ${OBJS_DRT_CONSTRAINT_LIB})
add_library(drt_lib_contact ${OBJS_DRT_CONTACT_LIB} ${OBJS_DRT_BEAMCONTACT_LIB})
add_library(drt_lib_beaminteraction ${OBJS_DRT_BEAMINTERACTION_LIB})
add_library(drt_lib_poromulti_art_coupl ${OBJS_DRT_POROMULTIPHASE_ART_COUPL})
add_library(drt_lib_xfem ${OBJS_DRT_XFEM})
add_library(drt_lib_io ${OBJS_DRT_IO_LIB})
#add_library(drt_lib_fluid ${OBJS_DRT_FLUID_LIB})
add_library(drt_lib_mortar ${OBJS_DRT_MORTAR_LIB})
add_library(drt_lib_nurbs ${OBJS_DRT_NURBS_DISCRET_LIB})
add_library(drt_lib_binstrategy ${OBJS_DRT_BINSTRATEGY_DISCRET_LIB})
add_library(drt_lib_particleengine ${OBJS_DRT_PARTICLE_ENGINE_LIB})
add_library(drt_lib_surfstress ${OBJS_DRT_SURFSTRESS})
add_library(drt_lib_utils ${OBJS_DRT_UTILS})
add_library(drt_lib_artnet ${OBJS_DRT_ARTNET_LIB})
add_library(drt_lib_redairways ${OBJS_DRT_RED_AIRWAYS_LIB})
add_library(drt_lib_cut ${OBJS_DRT_CUT})
add_library(drt_lib_statinv_ana ${OBJS_STATINV_ANALYSIS_LIB})
add_library(drt_lib_immersed ${OBJS_DRT_IMMERSED_LIB})
add_library(drt_lib_xfem_manag ${OBJS_DRT_XFEM_MANAG_LIB})
add_library(drt_lib_fiber ${OBJS_DRT_FIBER_LIB})
add_library(drt_lib_geometry_pair ${OBJS_DRT_GEOMETRY_PAIR})
add_library(alg_mor ${OBJS_DRT_MOR})


set(SUPPORT_LIBRARIES
  drt_lib_xfem_manag
  drt_lib
  drt_lib_artnet
  drt_lib_redairways
  drt_lib_mortar
  drt_lib_constraint
  drt_lib_contact
  drt_lib_beaminteraction
  drt_lib_poromulti_art_coupl
  # Libraries needed by elements
  drt_lib_surfstress
  drt_lib_xfem
  drt_lib_cut
  drt_lib_inpar
  drt_lib_io
  drt_lib_nurbs
  drt_lib_binstrategy
  drt_lib_particleengine
  drt_lib_statinv_ana
  drt_lib_immersed
  drt_lib_fiber
  drt_lib_geometry_pair
  )

# separate library as it does not play well with testing
add_library(drt_error ${OBJS_DRT_ERROR})

# filter framework
add_library(post_drt_common ${OBJS_POST_DRT_COMMON})
add_library(filter_common ${OBJS_FILTER_COMMON})
add_library(post_drt_ensight_lib ${OBJS_POST_DRT_ENSIGHT})
add_library(post_drt_vtk_lib ${OBJS_POST_DRT_VTK})

# linear algebra objects
#add_library(linalg SHARED ${OBJS_LINALG})
add_library(linalg ${OBJS_LINALG})
add_library(libsolver ${OBJS_SOLVER})

# continuation method
add_library(locacontinuation ${OBJS_LOCA_CONTINUATION})

# nonlinear solver objects
add_library(nlnsolver_nox ${OBJS_SOLVER_NONLIN_NOX})

string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
set(baciname baci-${BUILD_TYPE})

# the baci executable
add_executable (${baciname} ${OBJS_BACI_MAIN})

add_executable (post_processor EXCLUDE_FROM_ALL ${OBJS_POST_PROCESSOR_MAIN})
add_executable (post_drt_monitor EXCLUDE_FROM_ALL ${OBJS_POST_DRT_MONITOR})
add_executable (post_drt_s8convert EXCLUDE_FROM_ALL ${OBJS_POST_DRT_S8CONVERT})
add_executable (pre_exodus EXCLUDE_FROM_ALL ${OBJS_PRE_EXODUS})
add_executable (pre_locsys EXCLUDE_FROM_ALL ${OBJS_PRE_LOCSYS})

set(BACI_LIBRARIES
  ${MULTIFIELD_LIBRARIES}
  drt_adapter
  ${SINGLEFIELD_LIBRARIES}
  drt_nurbs
  drt_mortar
  drt_volmortar
  drt_constraint
  drt_cardiovascular0d
  drt_mor
  drt_lib_register
  ${ELEMENT_LIBRARIES}
  ${SUPPORT_LIBRARIES}
  drt_micromat
  drt_scatra_multi
  drt_mat
  drt_lib_utils
  locacontinuation
  nlnsolver_nox
  libsolver
  linalg
  drt_comm)

target_link_libraries (${baciname} ${BACI_LIBRARIES} ${LIBRARIES} drt_error)

set(FILTER_LIBRARIES
  post_drt_common
  drt_lib_register
  ${ELEMENT_LIBRARIES}
  ${SUPPORT_LIBRARIES}
  drt_mat
  drt_lib_utils
  libsolver
  linalg                  #linalg requires symbols from drt_lib which are in the SUPPORT_LIBRARIES
  filter_common
  drt_comm
  ${LIBRARIES}
  drt_error)

target_link_libraries (post_processor post_drt_ensight_lib post_drt_vtk_lib ${FILTER_LIBRARIES})
target_link_libraries (post_drt_monitor ${FILTER_LIBRARIES})
target_link_libraries (post_drt_s8convert ${FILTER_LIBRARIES})
target_link_libraries (pre_locsys ${FILTER_LIBRARIES})
target_link_libraries (pre_exodus ${FILTER_LIBRARIES})


ADD_CUSTOM_TARGET(baci DEPENDS ${baciname})
ADD_CUSTOM_TARGET(post DEPENDS post_processor post_drt_monitor)
ADD_CUSTOM_TARGET(post_drt_ensight DEPENDS post_processor)
ADD_CUSTOM_TARGET(post_drt_vtk DEPENDS post_processor)
ADD_CUSTOM_TARGET(framework DEPENDS pre_exodus post)
ADD_CUSTOM_TARGET(full DEPENDS ${baciname} post_processor post_drt_monitor post_drt_s8convert pre_exodus pre_locsys cut_test)

# various post-processing tools
ADD_CUSTOM_COMMAND(TARGET post_processor
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -D CMAKE_PROJECT_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR} -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/scripts/baci_create_post_scripts.cmake)

add_subdirectory(doc EXCLUDE_FROM_ALL)

###------------------------------------------------------------------ Tests

ENABLE_TESTING()
include(CTest)

###------------------------------------------------------------------ Unit Tests

add_executable(cut_test EXCLUDE_FROM_ALL ${OBJS_CUT_TEST})
target_link_libraries(cut_test drt_lib_cut ${QHULL_LIBRARIES} ${LAPACK_LIBRARIES}
    -Wl,--start-group ${FILTER_LIBRARIES} -Wl,--end-group
)
add_test(NAME run_cut_test COMMAND ${MPI_RUN} $<TARGET_FILE:cut_test>)

include (TestingFramework.cmake)
include (UnittestingFramework.cmake)
