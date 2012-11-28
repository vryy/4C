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
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/src/headers/compile_settings.h.in"
  "${PROJECT_BINARY_DIR}/src/headers/compile_settings.h"
  )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}/src/headers")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

# test for svn executable
FIND_PROGRAM(Subversion_SVN_EXECUTABLE svn
  DOC "subversion command line client")
MARK_AS_ADVANCED(Subversion_SVN_EXECUTABLE)

IF(Subversion_SVN_EXECUTABLE)

  # test if project source directory is a subversion directory
  # this test seems to be missing in the subversion package
  EXECUTE_PROCESS(COMMAND ${Subversion_SVN_EXECUTABLE} info ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE DUMMY_WC_INFO
    ERROR_VARIABLE Subversion_svn_info_error  
    RESULT_VARIABLE Subversion_svn_info_result
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  IF(${Subversion_svn_info_result} EQUAL 0)

    # get the package that does the real work
    FIND_PACKAGE(Subversion QUIET)
    IF(Subversion_FOUND)
      Subversion_WC_INFO(${CMAKE_HOME_DIRECTORY} baci)
      add_definitions("-DCHANGEDREVISION=${baci_WC_REVISION}")
      #add_definitions("-DCHANGEDDATE=${baci_WC_LAST_CHANGED_DATE}")
    ENDIF(Subversion_FOUND)
  ENDIF(${Subversion_svn_info_result} EQUAL 0)
ENDIF(Subversion_SVN_EXECUTABLE)

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
add_definitions("-DPARMETIS") # TODO remove PARMETIS flag

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

option (TRILINOS_DEV "Turn on advanced functionality from most recent Trilinos packages" OFF)
if (TRILINOS_DEV)
  add_definitions("-DTRILINOS_DEV")
endif (TRILINOS_DEV)

option (HAVE_TEKO "Turn on Teko support (needs Teko and Stratimikos to be compiled within Trilinos" OFF)
if (HAVE_TEKO)
  MESSAGE("II Turn on support for Teko")
  add_definitions("-DHAVE_TEKO")
endif (HAVE_TEKO)

option (HAVE_MueLu "Turn on MueLu support (needs MueLu to be compiled within Trilinos" OFF)
if (HAVE_MueLu)
  MESSAGE("II Turn on support for MueLu")
  add_definitions("-DHAVE_MueLu")
endif (HAVE_MueLu)

option (HAVE_Stratimikos "Turn on Stratimikos support (needs Stratimikos to be compiled within Trilinos" ON)
if (HAVE_Stratimikos)
  MESSAGE("II Turn on support for Stratimikos")
  add_definitions("-DHAVE_Stratimikos")
endif (HAVE_Stratimikos)

option (ENABLE_STACKTR "Enable printing of a stacktrace in case of a dserror" ON)
if (ENABLE_STACKTR)
  MESSAGE("II Print stacktrace in dserror")
  add_definitions("-DENABLE_STACKTR")
endif (ENABLE_STACKTR)

option (FS3I_INCACouplingTest "Nightly test for FS3I_INCA_coupling turned off" OFF)

#   additional features


  find_package(HDF5 REQUIRED)
  if (HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
    link_directories(${HDF5_LIBRARY_DIRS})
    set(LIBRARIES
      ${LIBRARIES}
      ${HDF5_LIBRARIES}
      )
  endif(HDF5_FOUND)


find_package(Qhull REQUIRED)
if (QHULL_FOUND)
  include_directories(${QHULL_INCLUDE_DIRS})
  link_directories(${QHULL_LIBRARY_DIRS})
  #MESSAGE("QHULL includes: ${QHULL_INCLUDE_DIRS}")
  #MESSAGE("QHULL library: ${QHULL_LIBRARY}")
  #MESSAGE("QHULL libraries: ${QHULL_LIBRARIES}")
  set(LIBRARIES
    ${LIBRARIES}
    ${QHULL_LIBRARIES}
    )
endif (QHULL_FOUND)

if (NOT QHULL_FOUND)
   MESSAGE(FATAL_ERROR "QHull not found")
endif ()

add_definitions("-DQHULL")

option (D_EXODUS "Pre_exodus (exodusII, netcdf)" OFF)
if (D_EXODUS)
  add_definitions("-DD_EXODUS")
endif(D_EXODUS)

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

# these are the fedora 8 with openmpi paths

if (NOT MPI_DIR)
  set (MPI_DIR "/usr")  # standard path for MPI in Fedora 8
endif (NOT MPI_DIR)

if (NOT INCLUDE_INSTALL_DIR)
  set ( INCLUDE_INSTALL_DIR "/lnm/lib/fc8/openmpi/include" )
endif (NOT INCLUDE_INSTALL_DIR)

if (NOT LIB_INSTALL_DIR)
  set ( LIB_INSTALL_DIR "/lnm/lib/fc8/openmpi/lib" )
endif (NOT LIB_INSTALL_DIR)

if (NOT Trilinos_DIR)
  set ( Trilinos_DIR "/lnm/lib/fc8/openmpi/trilinos" )
endif (NOT Trilinos_DIR)

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
IF (Trilinos_VERSION STREQUAL "10.1d") 
MESSAGE("   LNM DEFAULT")
ENDIF ()
MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
MESSAGE("End of Trilinos details\n")

FIND_PACKAGE( Boost COMPONENTS graph program_options REQUIRED)

if (Boost_FOUND)
  #message("Found Boost: ${Boost_INCLUDE_DIRS}")
  #message("Boost libraries: ${Boost_LIBRARIES}")
  #message("Boost libraries: ${Boost_LIBRARY_DIRS}")
  include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
  link_directories(${Boost_LIBRARY_DIRS}) 
  set(LIBRARIES
    ${Boost_LIBRARIES}
    ${LIBRARIES}
    )
  # Pass the version to baci
  add_definitions("-DBOOST_MAJOR_VERSION=${Boost_MAJOR_VERSION}")
  add_definitions("-DBOOST_MINOR_VERSION=${Boost_MINOR_VERSION}")
else()
 message(FATAL_ERROR ${Boost_ERROR_REASON})
endif()

# trilinos libraries
# Choose trilinos libraries explicitly.
if (Trilinos_FOUND)
  include_directories(SYSTEM ${Trilinos_INCLUDE_DIRS})
  include_directories(SYSTEM ${Trilinos_TPL_INCLUDE_DIRS})
  link_directories(${Trilinos_LIBRARY_DIRS})

  set(LIBRARIES
    ${Trilinos_LIBRARIES}
    ${Trilinos_TPL_LIBRARIES}
    ${LIBRARIES}
    )
    
# define of which library we are using - for compatibility checks
#if (Trilinos_VERSION STREQUAL "10.1d") 
#add_definitions("-DTRILINOS_10_1d")
#else ()
#endif ()
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

# FFTW for Jonas Biehler
find_package(FFTW)
if (FFTW_FOUND)
  include_directories(${FFTW_INCLUDE_DIRS})
  link_directories(${FFTW_LIBRARY_DIRS})
  set(LIBRARIES ${LIBRARIES} ${FFTW_LIBRARIES} )
  add_definitions("-DHAVE_FFTW")  # HAVE_FFTW preprocessor flag for BACI
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
set(OBJS_DRT_B2_LIB
  src/drt_beam2/beam2.cpp
  src/drt_beam2/beam2_evaluate.cpp
  src/drt_beam2/beam2_input.cpp
)
set(OBJS_DRT_LIB
  src/drt_fluid/drt_periodicbc.cpp
  src/drt_lib/drt_assemblestrategy.cpp
  src/drt_lib/drt_condition.cpp
  src/drt_lib/drt_conditiondefinition.cpp
  src/drt_lib/drt_container.cpp
  src/drt_lib/drt_discret.cpp
  src/drt_lib/drt_discret_xfem.cpp
  src/drt_lib/drt_discret_conditions.cpp
  src/drt_lib/drt_discret_evaluate.cpp
  src/drt_lib/drt_discret_fillcomplete.cpp
  src/drt_lib/drt_discret_partition.cpp
  src/drt_lib/drt_discret_utils.cpp
  src/drt_lib/drt_dofset.cpp
  src/drt_lib/drt_dofset_base.cpp
  src/drt_lib/drt_dofset_proxy.cpp
  src/drt_lib/drt_dofset_subproxy.cpp
  src/drt_lib/drt_dofset_transparent.cpp
  src/drt_lib/drt_dofset_transparent_independent.cpp
  src/drt_lib/drt_dofset_pbc.cpp
  src/drt_lib/drt_dofset_independent.cpp
  src/drt_lib/drt_dofset_independent_pbc.cpp
  src/drt_lib/drt_element.cpp
  src/drt_lib/drt_singletondestruction.cpp
  src/drt_lib/drt_function.cpp
  src/drt_lib/drt_globalproblem.cpp
  src/drt_lib/drt_locsys.cpp
  src/drt_lib/drt_node.cpp
  src/drt_lib/drt_parobject.cpp
  src/drt_lib/drt_parobjectfactory.cpp
  src/drt_lib/drt_timecurve.cpp
  src/drt_lib/drt_utils.cpp
  src/drt_lib/drt_utils_createdis.cpp
  src/drt_lib/drt_utils_factory.cpp
  src/drt_lib/drt_utils_nullspace.cpp
  src/drt_lib/drt_utils_parallel.cpp
)
set(OBJS_DRT_MAT_LIB
  src/drt_mat/aaagasser.cpp
  src/drt_mat/aaaneohooke.cpp
  src/drt_mat/aaaneohooke_stopro.cpp
  src/drt_mat/aaaraghavanvorp_damage.cpp
  src/drt_mat/aaa_mixedeffects.cpp
  src/drt_mat/anisotropic_balzani.cpp
  src/drt_mat/arrhenius_pv.cpp
  src/drt_mat/arrhenius_spec.cpp
  src/drt_mat/arrhenius_temp.cpp
  src/drt_mat/artwallremod.cpp
  src/drt_mat/biocell.cpp
  src/drt_mat/biofilm.cpp
  src/drt_mat/carreauyasuda.cpp
  src/drt_mat/charmm.cpp
  src/drt_mat/cnst_1d_art.cpp
  src/drt_mat/compogden.cpp
  src/drt_mat/constraintmixture.cpp
  src/drt_mat/contchainnetw.cpp
  src/drt_mat/elasthyper.cpp
  src/drt_mat/viscogenmax.cpp
  src/drt_mat/ferech_pv.cpp
  src/drt_mat/fluidporo.cpp
  src/drt_mat/fourieriso.cpp
  src/drt_mat/growth_ip.cpp
  src/drt_mat/holzapfelcardiovascular.cpp
  src/drt_mat/humphreycardiovascular.cpp
  src/drt_mat/ion.cpp
  src/drt_mat/itskov.cpp
  src/drt_mat/logneohooke.cpp
  src/drt_mat/material.cpp
  src/drt_mat/material_service.cpp
  src/drt_mat/matlist.cpp
  src/drt_mat/matpar_bundle.cpp
  src/drt_mat/matpar_material.cpp
  src/drt_mat/matpar_parameter.cpp
  src/drt_mat/micromaterial.cpp
  src/drt_mat/mixfrac.cpp
  src/drt_mat/modpowerlaw.cpp
  src/drt_mat/mooneyrivlin.cpp
  src/drt_mat/myocard.cpp
  src/drt_mat/neohooke.cpp
  src/drt_mat/newtonianfluid.cpp
  src/drt_mat/optimization_density.cpp
  src/drt_mat/permeablefluid.cpp
  src/drt_mat/plasticneohooke.cpp
  src/drt_mat/plastichyperelast.cpp
  src/drt_mat/plasticlinelast.cpp
  src/drt_mat/damage.cpp
  src/drt_mat/protein.cpp
  src/drt_mat/scatra_mat.cpp
  src/drt_mat/spring.cpp
  src/drt_mat/structporo.cpp
  src/drt_mat/stvenantkirchhoff.cpp
  src/drt_mat/sutherland.cpp
  src/drt_mat/thermostvenantkirchhoff.cpp
  src/drt_mat/thermoplasticlinelast.cpp
  src/drt_mat/robinson.cpp
  src/drt_mat/viscoanisotropic.cpp
  src/drt_mat/visconeohooke.cpp
  src/drt_mat/yeoh.cpp
  src/drt_mat/yoghurt.cpp
  src/drt_matelast/elast_coupanisoexpo.cpp
  src/drt_matelast/elast_coupanisoexpotwocoup.cpp
  src/drt_matelast/elast_coupanisoneohooke.cpp
  src/drt_matelast/elast_coupanisoneohooke_VarProp.cpp
  src/drt_matelast/elast_coupblatzko.cpp
  src/drt_matelast/elast_couplogneohooke.cpp
  src/drt_matelast/elast_coupneohooke.cpp
  src/drt_matelast/elast_coupvarga.cpp
  src/drt_matelast/elast_isoanisoexpo.cpp
  src/drt_matelast/elast_isocub.cpp
  src/drt_matelast/elast_isoexpopow.cpp
  src/drt_matelast/elast_coupmooneyrivlin.cpp
  src/drt_matelast/elast_isomooneyrivlin.cpp
  src/drt_matelast/elast_isoneohooke.cpp
  src/drt_matelast/elast_isoquad.cpp
  src/drt_matelast/elast_iso1pow.cpp
  src/drt_matelast/elast_iso2pow.cpp
  src/drt_matelast/elast_coup1pow.cpp
  src/drt_matelast/elast_coup2pow.cpp
  src/drt_matelast/elast_isovarga.cpp
  src/drt_matelast/elast_isoyeoh.cpp
  src/drt_matelast/elast_summand.cpp
  src/drt_matelast/elast_vologden.cpp
  src/drt_matelast/elast_volpenalty.cpp
  src/drt_matelast/elast_volsussmanbathe.cpp
  src/drt_matelast/elast_isovolaaagasser.cpp
  src/drt_matelast/elast_isovolHUdependentneohooke.cpp
  src/drt_patspec/patspec.cpp
  src/drt_mat/maxwell_0d_acinus.cpp
)
set(OBJS_DRT_STATMECH
  src/drt_statmech/statmech_manager.cpp
  src/drt_statmech/statmech_output.cpp
)
set(OBJS_DRT_ADAPT
  src/drt_adapter/adapter_algorithmbase.cpp
  src/drt_adapter/adapter_coupling.cpp
  src/drt_adapter/adapter_coupling_mortar.cpp
  src/drt_adapter/ad_fld_fluid.cpp
  src/drt_adapter/ad_fld_fluid_ale.cpp
  src/drt_adapter/ad_fld_base_algorithm.cpp
  src/drt_adapter/ad_fld_lung.cpp
  src/drt_adapter/ad_fld_moving_boundary.cpp
  src/drt_adapter/ad_fld_poro.cpp
  src/drt_adapter/ad_fld_fluid_fsi.cpp
  src/drt_adapter/ad_fld_fluid_fluid_fsi.cpp
  src/drt_adapter/ad_fld_fluid_xfem.cpp
  src/drt_adapter/ad_fld_xfluid_fsi.cpp
  src/drt_adapter/ad_opt.cpp
  src/drt_adapter/ad_opt_fluid_adjoint_base.cpp
  src/drt_adapter/ad_opt_fluid_algo.cpp
  src/drt_adapter/adapter_scatra_base_algorithm.cpp
  src/drt_adapter/adapter_scatra_fluid_ale_coupling_algo.cpp
  src/drt_adapter/adapter_scatra_fluid_coupling_algorithm.cpp
  src/drt_adapter/ad_str_structure.cpp
  src/drt_adapter/ad_str_bio.cpp
  src/drt_adapter/ad_str_constr_merged.cpp
  src/drt_adapter/ad_str_lung.cpp
  src/drt_adapter/ad_str_redairway.cpp
  src/drt_adapter/ad_str_timint_adaptive.cpp
  src/drt_adapter/ad_str_wrapper.cpp
  src/drt_adapter/ad_str_fsiwrapper.cpp
  src/drt_adapter/adapter_thermo.cpp
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
  src/drt_red_airways/airway_impl.cpp
  src/drt_red_airways/acinus_impl.cpp
  src/drt_red_airways/inter_acinar_dep_impl.cpp
  src/drt_red_airways/red_airway.cpp
  src/drt_red_airways/red_acinus.cpp
  src/drt_red_airways/red_inter_acinar_dep.cpp
  src/drt_red_airways/red_airway_input.cpp
  src/drt_red_airways/red_airway_resulttest.cpp
  src/drt_red_airways/redairway_tissue.cpp
)
set(OBJS_DRT_INPAR
  src/drt_inpar/drt_boolifyparameters.cpp
  src/drt_inpar/drt_validconditions.cpp
  src/drt_inpar/drt_validmaterials.cpp
  src/drt_inpar/drt_validparameters.cpp
  src/drt_inpar/inpar_problemtype.cpp
)
set(OBJS_DRT_ARTNET
  src/drt_art_net/art_net_dyn_drt.cpp
  src/drt_art_net/artnetexplicitintegration.cpp
)
set(OBJS_INTERSECTION
  src/drt_xfem/dof_distribution_switcher.cpp
  src/drt_xfem/dofkey.cpp
  src/drt_xfem/dof_management.cpp
  src/drt_xfem/dof_management_element.cpp
  src/drt_xfem/enrichment.cpp
  src/drt_xfem/enrichment_utils.cpp
  src/drt_xfem/field_enriched.cpp
  src/drt_xfem/physics.cpp
  src/drt_xfem/timeInt.cpp
  src/drt_xfem/timeInt_enr.cpp
  src/drt_xfem/timeInt_std_extrapolation.cpp
  src/drt_xfem/timeInt_std_SemiLagrange.cpp
  src/drt_xfem/xdofmapcreation_combust.cpp
  src/drt_xfem/xdofmapcreation_parallel_utils.cpp
  src/drt_xfem/xfem_edgestab.cpp
  src/drt_xfem/xfem_fluiddofset.cpp
  src/drt_xfem/xfem_fluidwizard.cpp
  src/drt_xfem/xfem_neumann.cpp
  src/drt_xfem/xfluidfluid_timeInt.cpp
  src/drt_xfem/xfluid_timeInt_base.cpp
  src/drt_xfem/xfluid_timeInt_std_SemiLagrange.cpp
  src/drt_xfem/xfluid_timeInt.cpp
  src/drt_geometry/geo_intersection.cpp
  src/drt_combust/combust_interface.cpp
  src/drt_combust/combust_flamefront.cpp
)
set(OBJS_POST_DRT_MONITOR
  src/post_monitor/post_drt_monitor.cpp
)
set(OBJS_DRT_ARTNET_LIB
  src/drt_art_net/art_junction.cpp
  src/drt_art_net/art_terminal_bc.cpp
  src/drt_art_net/art_write_gnuplot.cpp
  src/drt_art_net/artery.cpp
  src/drt_art_net/artery_evaluate.cpp
  src/drt_art_net/artery_input.cpp
  src/drt_art_net/artery_lin_exp.cpp
  src/drt_art_net/artery_resulttest.cpp
)
set(OBJS_DRT_COMBUST
  src/drt_combust/combust_algorithm.cpp
  src/drt_combust/combust_dyn.cpp
  src/drt_combust/combust_fluidimplicitintegration.cpp
  src/drt_combust/combust_fluidresulttest.cpp
  src/drt_combust/combust_reinitializer.cpp
  src/drt_combust/combust_reinitialization_pde.cpp
  src/drt_combust/combust_utils.cpp
)
set(OBJS_FILTER_COMMON
  src/filter_common/filter_evaluation.cpp
)
set(OBJS_DRT_BEAMCONTACT_LIB
  src/drt_beamcontact/beam3contact.cpp
  src/drt_beamcontact/beam3contact_manager.cpp
  src/drt_beamcontact/beam3contact_octtree.cpp
)
set(OBJS_DRT_STRU_MULTI
  src/drt_stru_multi/convert_mat.cpp
  src/drt_stru_multi/microstatic.cpp
  src/drt_stru_multi/microstatic_npsupport.cpp
  src/drt_stru_multi/microstatic_service.cpp
)
set(OBJS_SOH8_LIB
  src/drt_so3/inversedesign.cpp
  src/drt_so3/inversedesign_evaluate.cpp
  src/drt_so3/prestress.cpp
  src/drt_so3/so_disp.cpp
  src/drt_so3/so_disp_evaluate.cpp
  src/drt_so3/so_disp_input.cpp
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
  src/drt_so3/so_hex8_thermo_evaluate.cpp
  src/drt_so3/so_hex8p1j1.cpp
  src/drt_so3/so_hex8p1j1_evaluate.cpp
  src/drt_so3/so_hex8p1j1_input.cpp
  src/drt_so3/so_hex8fbar.cpp
  src/drt_so3/so_hex8fbar_evaluate.cpp
  src/drt_so3/so_hex8fbar_input.cpp
  src/drt_so3/so_line.cpp
  src/drt_so3/so_line_evaluate.cpp
  src/drt_so3/so_material.cpp
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
  src/drt_so3/so_tet10.cpp
  src/drt_so3/so_tet10_evaluate.cpp
  src/drt_so3/so_tet10_input.cpp
  src/drt_so3/so_tet4.cpp
  src/drt_so3/so_tet4_evaluate.cpp
  src/drt_so3/so_tet4_multiscale.cpp
  src/drt_so3/so_tet4_input.cpp
  src/drt_so3/so_weg6.cpp
  src/drt_so3/so_weg6_evaluate.cpp
  src/drt_so3/so_weg6_input.cpp
  src/drt_so3/so3_poro.cpp
  src/drt_so3/so3_poro_evaluate.cpp
  src/drt_so3/so3_poro_eletypes.cpp
  src/drt_so3/so3_scatra.cpp
  src/drt_so3/so3_scatra_evaluate.cpp
  src/drt_so3/so3_scatra_eletypes.cpp
  src/drt_so3/so3_thermo.cpp
  src/drt_so3/so3_thermo_evaluate.cpp
  src/drt_so3/so3_thermo_eletypes.cpp
)
set(OBJS_POST_DRT_ENSIGHT
  src/post_ensight/post_drt_ensight.cpp
  src/post_ensight/post_drt_ensight_single_field_writers.cpp
  src/post_ensight/post_drt_ensight_structure_stress.cpp
  src/post_ensight/post_drt_ensight_thermo_heatflux.cpp
  src/post_ensight/post_drt_ensight_writer.cpp
  src/post_ensight/post_drt_ensight_writer_nurbs.cpp
)
set(OBJS_POST_DRT_GENERIC
  src/post_generic/post_drt_generic.cpp
  src/post_generic/post_drt_generic_single_field_writers.cpp
  src/post_generic/post_drt_generic_structure_stress.cpp
  src/post_generic/post_drt_generic_writer.cpp
)
set(OBJS_DRT_ELCH
  src/drt_elch/elch_algorithm.cpp
  src/drt_elch/elch_dyn.cpp
  src/drt_elch/elch_moving_boundary_algorithm.cpp
)
set(OBJS_DRT_BELE3_LIB
  src/drt_bele3/bele2.cpp
  src/drt_bele3/bele2_evaluate.cpp
  src/drt_bele3/bele3_4.cpp
  src/drt_bele3/bele3_4_evaluate.cpp
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

  src/drt_fem_general/debug_nan.cpp
  src/drt_fem_general/drt_utils_boundary_integration.cpp
  src/drt_fem_general/drt_utils_fem_shapefunctions.cpp
  src/drt_fem_general/drt_utils_integration.cpp
  src/drt_fem_general/drt_utils_gausspoints.cpp
  src/drt_fem_general/drt_utils_local_connectivity_matrices.cpp
  src/drt_fem_general/largerotations.cpp
  src/drt_fluid/fluid_rotsym_periodicbc_utils.cpp
  src/drt_geometry/integrationcell.cpp
  src/drt_geometry/intersection_interfacepoint.cpp
  src/drt_geometry/intersection_math.cpp
  src/drt_geometry/intersection_service.cpp
  src/drt_geometry/position_array.cpp
  src/drt_geometry/searchtree.cpp
  src/drt_geometry/searchtree_geometry_service.cpp
  src/drt_geometry/searchtree_nearestobject.cpp
  src/drt_geometry/tetrahedradecomposition.cpp
  src/drt_io/io_gmsh.cpp
  src/drt_lib/drt_condition_selector.cpp
  src/drt_lib/drt_condition_utils.cpp
  src/drt_lib/drt_dirichletextractor.cpp
  src/drt_lib/drt_dserror.cpp
  src/drt_lib/drt_elementdefinition.cpp
  src/drt_lib/drt_elementtype.cpp
  src/drt_lib/drt_exporter.cpp
  src/drt_lib/drt_nodematchingoctree.cpp
  src/drt_lib/drt_parser.cpp
  src/drt_lib/drt_resulttest.cpp
  src/drt_lib/drt_utils_metis.cpp
  src/drt_lib/drt_utils_parmetis.cpp
  src/pss_full/pss_am.c
  src/pss_full/pss_table.c
  src/pss_full/pss_table_iter.c
)
set(OBJS_DRT_SCATRA
  src/drt_scatra/passive_scatra_algorithm.cpp
  src/drt_scatra/scatra_dyn.cpp
  src/drt_scatra/scatra_resulttest.cpp
  src/drt_scatra/scatra_timint_bdf2.cpp
  src/drt_scatra/scatra_timint_genalpha.cpp
  src/drt_scatra/scatra_timint_implicit.cpp
  src/drt_scatra/scatra_timint_implicit_service.cpp
  src/drt_scatra/scatra_timint_ost.cpp
  src/drt_scatra/scatra_timint_stat.cpp
  src/drt_scatra/scatra_utils.cpp
  src/drt_scatra/scatra_utils_clonestrategy.cpp
  src/drt_scatra/scatra_timint_tg.cpp
  src/drt_scatra/scatra_timint_reinitialization.cpp
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
  src/fortran/s8_fserv.f
  src/drt_s8/shell8/s8_mat_linel.c
  src/drt_s8/shell8/s8_mat_ogden.c
  src/drt_s8/shell8/s8_mat_ogden2.c
  src/drt_s8/shell8/s8_mat_ogden3.c
  src/drt_s8/shell8/s8_mattransform.c
  src/drt_s8/shell8/s8_tmat.c
)
set(OBJS_DRT_FSI
  src/drt_fsi/fs_monolithic.cpp
  src/drt_fsi/fsi_algorithm.cpp
  src/drt_fsi/fsi_constrmonolithic.cpp
  src/drt_fsi/fsi_constrmonolithic_fluidsplit.cpp
  src/drt_fsi/fsi_constrmonolithic_structuresplit.cpp
  src/drt_fsi/fsi_debugwriter.cpp
  src/drt_fsi/fsi_dirichletneumann.cpp
  src/drt_fsi/fsi_dirichletneumannslideale.cpp
  src/drt_fsi/fsi_dyn.cpp
  src/drt_fsi/fsi_fluid_ale.cpp
  src/drt_fsi/fsi_fluid_xfem.cpp
  src/drt_fsi/fsi_lagrangeprec.cpp
  src/drt_fsi/fsi_lungmonolithic.cpp
  src/drt_fsi/fsi_lungmonolithic_fluidsplit.cpp
  src/drt_fsi/fsi_lungmonolithic_structuresplit.cpp
  src/drt_fsi/fsi_matrixtransform.cpp
  src/drt_fsi/fsi_monolithic.cpp
  src/drt_fsi/fsi_monolithic_nonox.cpp
  src/drt_fsi/fsi_monolithic_linearsystem.cpp
  src/drt_fsi/fsi_monolithicfluidsplit.cpp
  src/drt_fsi/fsi_monolithicstructuresplit.cpp
  src/drt_fsi/fsi_fluidfluidmonolithic_structuresplit_nonox.cpp
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
  src/drt_fsi/fsi_lung_overlapprec.cpp
  src/drt_fsi/fsi_constr_overlapprec.cpp
  src/drt_fsi/fsi_overlapprec_fsiamg_analysis.cpp
  src/drt_fsi/fsi_overlapprec_fsiamg_smoothers.cpp
  src/drt_fsi/fsi_overlapprec_fsiamg_vcycles.cpp
  src/drt_fsi/fsi_partitioned.cpp
  src/drt_fsi/fsi_resulttest.cpp
  src/drt_fsi/fsi_statustest.cpp
  src/drt_fsi/fsi_structureale.cpp
  src/drt_fsi/fsi_utils.cpp
)
set(OBJS_DRT_FS3I
  src/drt_fs3i/fs3i_dyn.cpp
  src/drt_fs3i/fs3i_partitioned.cpp
  src/drt_fs3i/fs3i_partitioned_service.cpp
  src/drt_fs3i/fs3i_partitioned_1wc.cpp
  src/drt_fs3i/fs3i_partitioned_2wc.cpp
  src/drt_fs3i/biofilm_fsi.cpp
  src/drt_fs3i/aero_tfsi.cpp
  src/drt_fs3i/aero_tfsi_serv.cpp
)
set(OBJS_DRT_RED_AIRWAYS
  src/drt_red_airways/airwayimplicitintegration.cpp
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
set(OBJS_DRT_SCATRA_LIB
  src/drt_scatra/scatra_ele_boundary_impl.cpp
  src/drt_scatra/scatra_ele_impl.cpp
  src/drt_scatra/scatra_ele_impl_reinit.cpp
  src/drt_scatra/scatra_ele_impl_turbulence_service.cpp
  src/drt_scatra/scatra_ele_impl_utils.cpp
  src/drt_scatra/scatra_element.cpp
  src/drt_scatra/scatra_element_boundary_evaluate.cpp
  src/drt_scatra/scatra_element_evaluate.cpp
  src/drt_scatra/scatra_element_input.cpp
)
set(OBJS_DRT_MESHFREE_SCATRA_LIB
  src/drt_meshfree_discret/drt_meshfree_cell_utils.cpp
  src/drt_meshfree_discret/meshfree_scatra_cell.cpp
  src/drt_meshfree_discret/meshfree_scatra_cell_evaluate.cpp
  src/drt_meshfree_discret/meshfree_scatra_cell_boundary_evaluate.cpp
  src/drt_meshfree_discret/meshfree_scatra_impl.cpp
#  src/drt_meshfree_discret/meshfree_scatra_boundary_impl.cpp
)
set(OBJS_DRT_T2_LIB
  src/drt_truss2/truss2.cpp
  src/drt_truss2/truss2_evaluate.cpp
  src/drt_truss2/truss2_input.cpp
)
set(OBJS_DRT_COMBUST3_LIB
  src/drt_combust/combust_refinementcell.cpp
  src/drt_combust/combust3.cpp
  src/drt_combust/combust3_evaluate.cpp
  src/drt_combust/combust3_input.cpp
  src/drt_combust/combust3_line.cpp
  src/drt_combust/combust3_line_evaluate.cpp
  src/drt_combust/combust3_surface.cpp
  src/drt_combust/combust3_surface_evaluate.cpp
  src/drt_combust/combust3_sysmat.cpp
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
  src/drt_thermo/thr_contact.cpp
)
set(OBJS_POST_GID_TXT
  src/post_gid/gid_out_txt.c
)
set(OBJS_DRT_FLUID
  src/drt_fluid/drt_transfer_turb_inflow.cpp
  src/drt_fluid/dyn_smag.cpp
  src/drt_fluid/fluid_coupling_red_models.cpp
  src/drt_fluid/fluid_dyn_nln_drt.cpp
  src/drt_fluid/fluid_genalpha_integration.cpp
  src/drt_fluid/fluid_meshtying.cpp
  src/drt_fluid/fluid_discret_extractor.cpp
  src/drt_fluid/fluid_utils.cpp
  src/drt_fluid/fluid_utils_infnormscaling.cpp
  src/drt_fluid/fluid_utils_mapextractor.cpp
  src/drt_fluid/fluid_MHD_evaluate.cpp
  src/drt_fluid/fluid_windkessel_optimization.cpp
  src/drt_fluid/fluidimpedancecondition.cpp
  src/drt_fluid/fluid_volumetric_surfaceFlow_condition.cpp
  src/drt_fluid/fluidimplicitintegration.cpp
  src/drt_fluid/fluidresulttest.cpp
  src/drt_fluid/scale_sep_gmo.cpp
  src/drt_fluid/time_integration_scheme.cpp
  src/drt_fluid/turbulent_flow_algorithm.cpp
  src/drt_fluid/turbulence_statistic_manager.cpp
  src/drt_fluid/turbulence_statistics_bfs.cpp
  src/drt_fluid/turbulence_statistics_oracles.cpp
  src/drt_fluid/turbulence_statistics_ccy.cpp
  src/drt_fluid/turbulence_statistics_cha.cpp
  src/drt_fluid/turbulence_statistics_bcf.cpp
  src/drt_fluid/turbulence_statistics_ldc.cpp
  src/drt_fluid/turbulence_statistics_mean_general.cpp
  src/drt_fluid/turbulence_statistics_sqc.cpp
  src/drt_fluid/xfluid.cpp
  src/drt_fluid/xfluidfluid.cpp
  src/drt_fluid/xfluidfluidresulttest.cpp
  src/drt_fluid/fluid_timint.cpp
)
set(OBJS_DRT_SR_LIB
  src/drt_smoothrod/smoothrod.cpp
  src/drt_smoothrod/smoothrod_evaluate.cpp
  src/drt_smoothrod/smoothrod_input.cpp
)
set(OBJS_DRT_TO2_LIB
  src/drt_torsion2/torsion2.cpp
  src/drt_torsion2/torsion2_evaluate.cpp
  src/drt_torsion2/torsion2_input.cpp
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
set(OBJS_POST_DRT_GID
  src/post_gid/post_drt_gid.cpp
)
set(OBJS_DRT_NURBS_DISCRET_LIB
  src/drt_fem_general/drt_utils_bspline.cpp
  src/drt_fem_general/drt_utils_nurbs_shapefunctions.cpp
  src/drt_nurbs_discret/drt_control_point.cpp
  src/drt_nurbs_discret/drt_knotvector.cpp
  src/drt_nurbs_discret/drt_nurbs_discret.cpp
)
set(OBJS_DRT_MESHFREE_DISCRET_LIB
  src/drt_fem_general/drt_utils_maxent_basisfunctions.cpp
  src/drt_meshfree_discret/drt_meshfree_node.cpp
  src/drt_meshfree_discret/drt_meshfree_bin.cpp
  src/drt_meshfree_discret/drt_meshfree_cell.cpp
  src/drt_meshfree_discret/drt_meshfree_cell_utils.cpp
  src/drt_meshfree_discret/drt_meshfree_discret.cpp
  src/drt_meshfree_discret/drt_meshfree_discret_evaluate.cpp
  src/drt_meshfree_discret/drt_meshfree_discret_fillcomplete.cpp
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
set(OBJS_POST_DRT_COMMON
  src/post_drt_common/post_drt_common.cpp
)
set(OBJS_INV_ANALYSIS
  src/drt_inv_analysis/gen_inv_analysis.cpp
  src/drt_inv_analysis/inv_analysis.cpp
  src/drt_structure/str_invanalysis.cpp
)  
set(OBJS_MLMC
  src/drt_mlmc/mlmc.cpp
  src/drt_structure/str_mlmc.cpp 
  src/drt_mlmc/gen_randomfield.cpp
)
set(OBJS_DRT_MAT
  src/drt_mat/micromaterial_evaluate.cpp
  src/drt_mat/micromaterialgp_static.cpp
)
set(OBJS_BACI_MAIN
  src/global_full/global_init_control.cpp
  src/global_full/global_inp_control.cpp
  src/global_full/global_cal_control.cpp
  src/global_full/global_control.cpp
  src/global_full/baci.cpp
)
set(OBJS_DRT_B3II_LIB
  src/drt_beam3ii/beam3ii.cpp
  src/drt_beam3ii/beam3ii_evaluate.cpp
  src/drt_beam3ii/beam3ii_input.cpp
)

set(OBJS_DRT_B3CL_LIB
  src/drt_beam3cl/beam3cl.cpp
  src/drt_beam3cl/beam3cl_evaluate.cpp
  src/drt_beam3cl/beam3cl_input.cpp
)
set(OBJS_DRT_B3EB_LIB
  src/drt_beam3eb/beam3eb.cpp
  src/drt_beam3eb/beam3eb_evaluate.cpp
  src/drt_beam3eb/beam3eb_input.cpp
)
set(OBJS_DRT_B3EBTOR_LIB
  src/drt_beam3ebtor/beam3ebtor.cpp
  src/drt_beam3ebtor/beam3ebtor_evaluate.cpp
  src/drt_beam3ebtor/beam3ebtor_input.cpp
)
set(OBJS_PRE_CURVE
  src/pre_curve/pre_curve.cpp
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
  src/drt_mortar/mortar_projector.cpp
  src/drt_mortar/mortar_strategy_base.cpp
)
set(OBJS_DRT_MORTAR_LIB
  src/drt_mortar/mortar_element.cpp
  src/drt_mortar/mortar_element_integrator.cpp
  src/drt_mortar/mortar_element_shapefct.cpp
  src/drt_mortar/mortar_node.cpp
)
set(OBJS_DRT_STRUCTURE
  src/drt_structure/strtimada.cpp
  src/drt_structure/strtimada_create.cpp
  src/drt_structure/strtimada_zienxie.cpp
  src/drt_structure/strtimint.cpp
  src/drt_structure/strtimint_create.cpp
  src/drt_structure/strtimint_ab2.cpp
  src/drt_structure/strtimint_centrdiff.cpp
  src/drt_particle/particle_timint_centrdiff.cpp
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
  src/drt_structure/stru_dyn_nln_drt.cpp
  src/drt_structure/stru_resulttest.cpp
  src/drt_structure/strtimint_statmech.cpp
)
set(OBJS_POST_DRT_S8CONVERT
  src/post_s8convert/post_drt_s8convert.cpp
)
set(OBJS_DRT_PARTICLE
	src/drt_particle/particle_algorithm.cpp
	src/drt_particle/particle_dyn.cpp
	src/drt_particle/cavitation_algorithm.cpp
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
)
set(OBJS_SOLVER
  src/linalg/linalg_solver.cpp
  src/solver/amg/aggregation_method.cpp
  src/solver/amg/aggregation_method_ml.cpp
  src/solver/amg/aggregation_method_uncoupled.cpp
  src/solver/amg/smoother.cpp
  src/solver/amg/transfer_operator.cpp
  src/solver/amg/transfer_operator_pgamg.cpp
  src/solver/amg/transfer_operator_pgamg2.cpp
  src/solver/amg/transfer_operator_saamg.cpp
  src/solver/amg/transfer_operator_tentative.cpp
  src/solver/amg/amgpreconditioner.cpp
  src/solver/amg/braesssarazin_smoother.cpp
  src/solver/saddlepointpreconditioner.cpp
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
  src/solver/solver_stratimikossolver.cpp
  src/solver/teko/teko_baciepetraoperatorwrapper.cpp
  src/solver/muelu/muelu_ContactAFilterFactory.cpp
  src/solver/muelu/muelu_ContactTransferFactory.cpp
  src/solver/muelu/muelu_ContactInfoFactory.cpp
  src/solver/muelu/muelu_ContactASlaveDofFilterFactory.cpp 
  src/solver/muelu/muelu_ContactSPAggregationFactory.cpp
  src/solver/muelu/MueLu_MyTrilinosSmoother.cpp
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
)
set(OBJS_DRT_IO_LIB
  src/drt_io/io.cpp
  src/drt_io/io_control.cpp
  src/drt_io/io_hdf.cpp
  src/drt_io/io_pstream.cpp
  src/drt_lib/drt_elementreader.cpp
  src/drt_lib/drt_inputreader.cpp
  src/drt_lib/drt_linedefinition.cpp
  src/drt_lib/drt_materialdefinition.cpp
  src/drt_lib/drt_nodereader.cpp
)
set(OBJS_FORTRAN
  src/fortran/hex_fserv.f
)
set(OBJS_DRT_B2R_LIB
  src/drt_beam2r/beam2r.cpp
  src/drt_beam2r/beam2r_evaluate.cpp
  src/drt_beam2r/beam2r_input.cpp
)
set(OBJS_DRT_POTENTIAL
  src/drt_potential/drt_potential.cpp
  src/drt_potential/drt_potential_container.cpp
  src/drt_potential/drt_potential_dofset.cpp
  src/drt_potential/drt_potential_manager.cpp
  src/drt_potential/drt_potential_surface.cpp
  src/drt_potential/drt_potential_volume.cpp
)
set(OBJS_DRT_NURBS_DISCRET
  src/drt_nurbs_discret/drt_apply_nurbs_initial_condition.cpp
  src/drt_nurbs_discret/drt_nurbs_utils.cpp
)
#set(OBJS_DRT_MESHFREE_DISCRET
#  src/drt_meshfree_discret/drt_meshfree_utils.cpp
#)
set(OBJS_DRT_B3_LIB
  src/drt_beam3/beam3.cpp
  src/drt_beam3/beam3_evaluate.cpp
  src/drt_beam3/beam3_input.cpp
)
set(OBJS_DRT_FLUID_ELE_LIB
  src/drt_fluid_ele/fluid_ele.cpp
  src/drt_fluid_ele/fluid_ele_boundary_calc.cpp
  src/drt_fluid_ele/fluid_ele_boundary.cpp
  src/drt_fluid_ele/fluid_ele_boundary_evaluate.cpp
  src/drt_fluid_ele/fluid_ele_calc.cpp
  src/drt_fluid_ele/fluid_ele_calc_ale_service.cpp
  src/drt_fluid_ele/fluid_ele_calc_general_service.cpp
  src/drt_fluid_ele/fluid_ele_calc_loma.cpp
  src/drt_fluid_ele/fluid_ele_intfaces.cpp
  src/drt_fluid_ele/fluid_ele_calc_intfaces_stab.cpp
  src/drt_fluid_ele/fluid_ele_intfaces_calc.cpp
  src/drt_fluid_ele/fluid_ele_calc_loma_service.cpp
  src/drt_fluid_ele/fluid_ele_calc_poro.cpp
  src/drt_fluid_ele/fluid_ele_calc_stabilization.cpp
  src/drt_fluid_ele/fluid_ele_calc_std.cpp
  src/drt_fluid_ele/fluid_ele_calc_turbulence_service.cpp
  src/drt_fluid_ele/fluid_ele_calc_weak_dbc.cpp
  src/drt_fluid_ele/fluid_ele_calc_xfem.cpp
  src/drt_fluid_ele/fluid_ele_calc_xfem_coupling.cpp
  src/drt_fluid_ele/fluid_ele_calc_xfem_coupling_impl.cpp
  src/drt_fluid_ele/fluid_ele_evaluate.cpp
  src/drt_fluid_ele/fluid_ele_factory.cpp
  src/drt_fluid_ele/fluid_ele_input.cpp
  src/drt_fluid_ele/fluid_ele_parameter.cpp
  src/drt_fluid_ele/fluid_genalpha_resVMM.cpp
  src/drt_fluid/time_integration_element.cpp
  src/drt_opti/topopt_fluidAdjoint3_boundary.cpp
  src/drt_opti/topopt_fluidAdjoint3_impl.cpp
  src/drt_opti/topopt_fluidAdjoint3_impl_parameter.cpp
)
set(OBJS_DRT_T3_LIB
  src/drt_truss3/truss3.cpp
  src/drt_truss3/truss3_evaluate.cpp
  src/drt_truss3/truss3_input.cpp
)
set(OBJS_DRT_ALEALG
  src/drt_ale/ale.cpp
  src/drt_ale/ale_dyn.cpp
  src/drt_ale/ale_laplace.cpp
  src/drt_ale/ale_lin.cpp
  src/drt_ale/ale_resulttest.cpp
  src/drt_ale/ale_springs.cpp
  src/drt_ale/ale_springs_fixed_ref.cpp
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
)
set(OBJS_DRT_CONTACT
  src/drt_contact/contact_abstract_strategy.cpp
  src/drt_contact/contact_coupling2d.cpp
  src/drt_contact/contact_coupling3d.cpp
  src/drt_contact/contact_integrator.cpp
  src/drt_contact/contact_interface.cpp
  src/drt_contact/contact_interface_tools.cpp
  src/drt_contact/contact_lagrange_strategy.cpp
  src/drt_contact/contact_manager.cpp
  src/drt_contact/contact_penalty_strategy.cpp
  src/drt_contact/meshtying_abstract_strategy.cpp
  src/drt_contact/meshtying_lagrange_strategy.cpp
  src/drt_contact/meshtying_manager.cpp
  src/drt_contact/meshtying_penalty_strategy.cpp
  src/drt_contact/selfcontact_binarytree.cpp
)
set(OBJS_DRT_CUT
  src/drt_cut/base_vol.cpp
  src/drt_cut/cut_boundarycell.cpp
  src/drt_cut/cut_boundingbox.cpp
  src/drt_cut/cut_coloredgraph.cpp
  src/drt_cut/cut_creator.cpp
  src/drt_cut/cut_cycle.cpp
  src/drt_cut/direct_divergence.cpp
  src/drt_cut/cut_edge.cpp
  src/drt_cut/cut_element.cpp
  src/drt_cut/cut_elementhandle.cpp
  src/drt_cut/cut_facet.cpp
  src/drt_cut/cut_facetgraph.cpp
  src/drt_cut/cut_find_cycles.cpp
  src/drt_cut/cut_integrationcell.cpp
  src/drt_cut/cut_integrationcellcreator.cpp
  src/drt_cut/cut_kernel.cpp
  src/drt_cut/cut_levelsetintersection.cpp
  src/drt_cut/cut_levelsetside.cpp
  src/drt_cut/cut_line.cpp
  src/drt_cut/cut_mesh.cpp
  src/drt_cut/cut_meshhandle.cpp
  src/drt_cut/cut_meshintersection.cpp
  src/drt_cut/cut_node.cpp
  src/drt_cut/cut_parallel.cpp
  src/drt_cut/cut_point.cpp
  src/drt_cut/cut_pointgraph.cpp
  src/drt_cut/cut_pointpool.cpp
  src/drt_cut/cut_side.cpp
  src/drt_cut/cut_sidehandle.cpp
  src/drt_cut/cut_tetmesh.cpp
  src/drt_cut/cut_tetmeshintersection.cpp
  src/drt_cut/cut_triangulateFacet.cpp
  src/drt_cut/cut_volumecell.cpp
  src/drt_cut/cut_volumecellgenerator.cpp
  src/drt_cut/line_integration.cpp
  src/drt_cut/facet_integration.cpp
  src/drt_cut/boundarycell_integration.cpp
  src/drt_cut/volume_integration.cpp
  src/drt_cut/least_squares.cpp
)
set(OBJS_CUT_TEST
  tests/cut_test/cut_test_double.cpp
  tests/cut_test/cut_test_facets.cpp
  tests/cut_test/cut_test_graph.cpp
  tests/cut_test/cut_test_generator.cpp
  tests/cut_test/cut_test_geometry.cpp
  tests/cut_test/cut_test_intersection.cpp
  tests/cut_test/cut_test_levelset.cpp
  tests/cut_test/cut_test_loader.cpp
  tests/cut_test/cut_test_main.cpp
  tests/cut_test/cut_test_simple.cpp
  tests/cut_test/cut_test_triangulateFacet.cpp
  tests/cut_test/cut_test_units.cpp
  tests/cut_test/cut_test_utils.cpp
  tests/cut_test/cut_test_fluidfluid.cpp
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
  src/drt_fem_general/drt_utils_gausspoints.cpp
  src/drt_fem_general/drt_utils_fem_shapefunctions.cpp
  src/drt_fem_general/drt_utils_integration.cpp
)
set(OBJS_GRAPH_TEST
  tests/graph_test/graph_test_main.cpp
  tests/graph_test/graph_test_read.cpp
  tests/graph_test/graph_test_gnuplot.cpp
  tests/graph_test/graph_test_output.cpp
  tests/graph_test/graph_test_find_cycles.cpp
)
set(OBJS_DRT_STRU_ALE
  src/drt_stru_ale/stru_ale_algorithm.cpp  
  src/drt_stru_ale/stru_ale_dyn.cpp
)  
set(OBJS_DRT_POROELAST
  src/drt_poroelast/poro_dyn.cpp
  src/drt_poroelast/poro_base.cpp
  src/drt_poroelast/poroelast_monolithic.cpp
  src/drt_poroelast/poroelast_utils.cpp
  src/drt_poroelast/poro_partitioned.cpp
  src/drt_poroelast/poro_monolithicsplit.cpp
  src/drt_poroelast/poro_monolithicstructuresplit.cpp
  src/drt_poroelast/poro_monolithicfluidsplit.cpp
  src/drt_poroelast/poro_scatra.cpp
)

set(OBJS_DRT_SSI
src/drt_ssi/ssi_dyn.cpp
src/drt_ssi/ssi_base.cpp
src/drt_ssi/ssi_partitioned.cpp
src/drt_ssi/ssi_partitioned_1wc.cpp
src/drt_ssi/ssi_partitioned_2wc.cpp
)
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

option (D_ARTNET      "1D Arterial network elements" OFF)
if (D_ARTNET)
  add_definitions("-DD_ARTNET")
endif(D_ARTNET)

option (D_ALE_BFLOW   "ALE deformation included in flow calculation" OFF)
if (D_ALE_BFLOW)
  add_definitions("-DD_ALE_BFLOW")
endif(D_ALE_BFLOW)

option (D_RED_AIRWAYS "reduced dimensional airways for calculating flow in lungs" OFF)
if (D_RED_AIRWAYS)
  add_definitions("-DD_RED_AIRWAYS")
endif(D_RED_AIRWAYS)

###------------------------------------------------------------------ Element flags

#   other flags
option (HAVENOT_SUPERLU "if no super lu is available, turn it on" OFF)
if (HAVENOT_SUPERLU)
  add_definitions("-DHAVENOT_SUPERLU")
endif(HAVENOT_SUPERLU)

###------------------------------------------------------------------ Binaries

# hide libraries

if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/src)
endif (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

# Define the libraries. The link order is important!

# Compile element libraries first. There are files that take time to compile.
add_library(ele_fluid ${OBJS_DRT_FLUID_ELE_LIB})
add_library(ele_ale ${OBJS_DRT_ALE_LIB})
add_library(ele_b2 ${OBJS_DRT_B2_LIB})
add_library(ele_b2r ${OBJS_DRT_B2R_LIB})
add_library(ele_b3 ${OBJS_DRT_B3_LIB})
add_library(ele_b3ii ${OBJS_DRT_B3II_LIB})
add_library(ele_b3eb ${OBJS_DRT_B3EB_LIB})
add_library(ele_b3ebtor ${OBJS_DRT_B3EBTOR_LIB})
add_library(ele_b3cl ${OBJS_DRT_B3CL_LIB})
add_library(ele_bele3 ${OBJS_DRT_BELE3_LIB})
add_library(ele_combust3 ${OBJS_DRT_COMBUST3_LIB})
add_library(ele_opti ${OBJS_DRT_TOPOPT_LIB})
add_library(ele_s8 ${OBJS_DRT_S8_LIB})
add_library(ele_scatra ${OBJS_DRT_SCATRA_LIB})
add_library(ele_meshfree_scatra ${OBJS_DRT_MESHFREE_SCATRA_LIB})
add_library(ele_so ${OBJS_SOH8_LIB})
add_library(ele_t2 ${OBJS_DRT_T2_LIB})
add_library(ele_t3 ${OBJS_DRT_T3_LIB})
add_library(ele_thermo ${OBJS_DRT_THERMO_LIB})
add_library(ele_to2 ${OBJS_DRT_TO2_LIB})
add_library(ele_to3 ${OBJS_DRT_TO3_LIB})
add_library(ele_w1 ${OBJS_DRT_W1_LIB})


set(ELEMENT_LIBRARIES 
  ele_ale
  ele_b2
  ele_b2r
  ele_b3
  ele_b3ii
  ele_b3eb
  ele_b3ebtor
  ele_b3cl
  ele_bele3
  ele_combust3
  ele_fluid
  ele_opti
  ele_s8
  ele_scatra
  ele_meshfree_scatra
  ele_so
  ele_t2
  ele_t3
  ele_thermo
  ele_to2
  ele_to3
  ele_w1
  )

# algorithm libraries that depend on adapters
add_library(alg_elch ${OBJS_DRT_ELCH})
add_library(alg_fsi ${OBJS_DRT_FSI})
add_library(alg_loma ${OBJS_DRT_LOMA})
add_library(alg_tsi ${OBJS_DRT_TSI})
add_library(alg_stru_ale ${OBJS_DRT_STRU_ALE})
add_library(alg_poroelast ${OBJS_DRT_POROELAST})
add_library(alg_fs3i ${OBJS_DRT_FS3I})
add_library(alg_topopt ${OBJS_DRT_TOPOPT})
add_library(alg_ssi ${OBJS_DRT_SSI})
add_library(alg_particle ${OBJS_DRT_PARTICLE})

set(MULTIFIELD_LIBRARIES
  alg_elch
  alg_fsi
  alg_loma
  alg_tsi
  alg_stru_ale
  alg_poroelast
  alg_fs3i
  alg_topopt
  alg_ssi
  alg_particle
  )

# the adapter layer
add_library(drt_adapter ${OBJS_DRT_ADAPT})

# algorithm libraries the adapters depend upon
add_library(alg_ale ${OBJS_DRT_ALEALG})
add_library(alg_artnet ${OBJS_DRT_ARTNET})
add_library(alg_combust ${OBJS_DRT_COMBUST})
#add_library(alg_contact SHARED ${OBJS_DRT_CONTACT})
add_library(alg_contact ${OBJS_DRT_CONTACT})
add_library(alg_fluid ${OBJS_DRT_FLUID})
add_library(alg_invanalysis ${OBJS_INV_ANALYSIS})
add_library(alg_mlmc ${OBJS_MLMC})
add_library(alg_redairways ${OBJS_DRT_RED_AIRWAYS})
add_library(alg_scatra ${OBJS_DRT_SCATRA})
add_library(alg_statmech ${OBJS_DRT_STATMECH})
add_library(alg_structure ${OBJS_DRT_STRUCTURE})
add_library(alg_thermo ${OBJS_DRT_THERMO})
add_library(alg_fluid_adjoint ${OBJS_DRT_TOPOPT_ADJOINT})
add_library(alg_opti ${OBJS_OPTI})

set(SINGLEFIELD_LIBRARIES
  alg_ale
  alg_artnet
  alg_combust
  alg_fluid
  alg_redairways
  alg_scatra
  alg_structure
  alg_statmech
  alg_invanalysis
  alg_mlmc
  alg_thermo
  alg_contact
  alg_fluid_adjoint
  alg_opti
  )

# general framework
add_library(drt_nurbs ${OBJS_DRT_NURBS_DISCRET})
#add_library(drt_meshfree ${OBJS_DRT_MESHFREE_DISCRET})
add_library(drt_constraint ${OBJS_DRT_CONSTRAINT})
#add_library(drt_mortar SHARED ${OBJS_DRT_MORTAR})
add_library(drt_mortar ${OBJS_DRT_MORTAR})

# material definitions
add_library(drt_mat ${OBJS_DRT_MAT_LIB})
add_library(drt_micromat ${OBJS_DRT_MAT} ${OBJS_DRT_STRU_MULTI})

# special library that demands all parobject types
add_library(drt_lib_register ${OBJS_DRT_REGISTER})

# libraries that are used by baci and by filters
add_library(drt_comm ${OBJS_DRT_COMM})
add_library(drt_lib ${OBJS_DRT_LIB})
add_library(drt_lib_inpar ${OBJS_DRT_INPAR})
add_library(drt_lib_constraint ${OBJS_DRT_CONSTRAINT_LIB})
add_library(drt_lib_contact ${OBJS_DRT_CONTACT_LIB} ${OBJS_DRT_BEAMCONTACT_LIB})
add_library(drt_lib_intersection ${OBJS_INTERSECTION})
add_library(drt_lib_io ${OBJS_DRT_IO_LIB})
#add_library(drt_lib_fluid ${OBJS_DRT_FLUID_LIB})
add_library(drt_lib_mortar ${OBJS_DRT_MORTAR_LIB})
add_library(drt_lib_nurbs ${OBJS_DRT_NURBS_DISCRET_LIB})
add_library(drt_lib_meshfree ${OBJS_DRT_MESHFREE_DISCRET_LIB})
add_library(drt_lib_potential ${OBJS_DRT_POTENTIAL})
add_library(drt_lib_surfstress ${OBJS_DRT_SURFSTRESS})
add_library(drt_lib_utils ${OBJS_DRT_UTILS})
add_library(drt_lib_artnet ${OBJS_DRT_ARTNET_LIB})
add_library(drt_lib_redairways ${OBJS_DRT_RED_AIRWAYS_LIB})
add_library(drt_lib_cut ${OBJS_DRT_CUT})

set(SUPPORT_LIBRARIES
  drt_lib
  drt_lib_artnet 
  drt_lib_redairways
  drt_lib_mortar
  drt_lib_constraint
  drt_lib_contact
  # Libraries needed by elements
  drt_lib_surfstress
  drt_lib_potential
  drt_lib_intersection
  drt_lib_cut
  drt_lib_inpar
  drt_lib_io
  drt_lib_nurbs 
  drt_lib_meshfree
  )

# filter framework
add_library(post_drt_common ${OBJS_POST_DRT_COMMON})
add_library(filter_common ${OBJS_FILTER_COMMON})

# linear algebra objects
#add_library(linalg SHARED ${OBJS_LINALG})
add_library(linalg ${OBJS_LINALG})
add_library(libsolver ${OBJS_SOLVER})

string(TOLOWER ${CMAKE_BUILD_TYPE} BUILD_TYPE)
set(baciname baci-${BUILD_TYPE})
#set(baciname baci)

# the baci executable
add_executable (${baciname} ${OBJS_BACI_MAIN})

add_executable (post_drt_gid     ${OBJS_POST_DRT_GID} ${OBJS_POST_GID_TXT})
add_executable (post_drt_ensight ${OBJS_POST_DRT_ENSIGHT})
add_executable (post_drt_generic ${OBJS_POST_DRT_GENERIC})
add_executable (post_drt_monitor ${OBJS_POST_DRT_MONITOR})
add_executable (post_drt_s8convert ${OBJS_POST_DRT_S8CONVERT})

add_executable (pre_curve ${OBJS_PRE_CURVE})

target_link_libraries (${baciname} 
  ${MULTIFIELD_LIBRARIES} 
  drt_adapter 
  ${SINGLEFIELD_LIBRARIES}
  drt_nurbs
#  drt_meshfree
  drt_mortar 
  drt_constraint
  drt_lib_register
  ${ELEMENT_LIBRARIES}
  ${SUPPORT_LIBRARIES}
  drt_micromat drt_mat
  drt_lib_utils
  libsolver
  linalg
  drt_comm
  ${LIBRARIES})

set(FILTER_LIBRARIES 
  post_drt_common
  drt_lib_register
  ${ELEMENT_LIBRARIES}
  ${SUPPORT_LIBRARIES}
  drt_mat
  drt_lib_utils
  libsolver
  linalg
  filter_common
  drt_comm
  ${LIBRARIES})

target_link_libraries (post_drt_gid ${FILTER_LIBRARIES})
target_link_libraries (post_drt_ensight ${FILTER_LIBRARIES})
target_link_libraries (post_drt_generic ${FILTER_LIBRARIES})
target_link_libraries (post_drt_monitor ${FILTER_LIBRARIES})
target_link_libraries (post_drt_s8convert ${FILTER_LIBRARIES})
target_link_libraries (pre_curve ${FILTER_LIBRARIES})

option (D_EXODUS "exodus filter" OFF)
if (D_EXODUS)
  add_definitions("-DD_EXODUS")
  add_executable (pre_exodus ${OBJS_PRE_EXODUS})
  target_link_libraries (pre_exodus ${FILTER_LIBRARIES})
endif (D_EXODUS)

###------------------------------------------------------------------ Tests

ENABLE_TESTING()
include(CTest)

###------------------------------------------------------------------ Unit Tests

if(NOT TRILINOS_DEV)
  if (NOT ACML_FOUND )  ## cut_test not available on kaiser
    add_executable(cut_test ${OBJS_CUT_TEST})
    target_link_libraries(cut_test drt_lib_cut ${QHULL_LIBRARIES} shards epetra teuchos ${LAPACK_LIBRARIES}
    drt_lib_io
    ${LIBRARIES}
    )
    add_test(NAME run_cut_test COMMAND ${MPI_DIR}/bin/mpirun $<TARGET_FILE:cut_test>)
  endif(NOT ACML_FOUND)
else()
# linker command for Q3/2012 Trilinos or newer
  if (NOT ACML_FOUND )  ## cut_test not available on kaiser
    add_executable(cut_test ${OBJS_CUT_TEST})
    target_link_libraries(cut_test drt_lib_cut ${QHULL_LIBRARIES} shards epetra teuchoscore ${LAPACK_LIBRARIES}
    drt_lib_io
    ${LIBRARIES}
    )
    add_test(NAME run_cut_test COMMAND ${MPI_DIR}/bin/mpirun $<TARGET_FILE:cut_test>)
  endif(NOT ACML_FOUND)
endif(NOT TRILINOS_DEV)

#add_executable(graph_test ${OBJS_GRAPH_TEST})
#target_link_libraries(graph_test ${Boost_LIBRARIES})

###------------------------------------------------------------------ Test definitions

macro (baci_test arg nproc restart)
  add_test(NAME ${arg}-p${nproc}
    COMMAND ${MPI_DIR}/bin/mpirun -np ${nproc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${arg}.dat xxx)

  if (${restart})
    add_test(NAME ${arg}-p${nproc}-restart
      COMMAND ${MPI_DIR}/bin/mpirun -np ${nproc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${arg}.dat xxx restart=${restart})
  endif (${restart})

  set_tests_properties ( ${arg}-p${nproc} PROPERTIES TIMEOUT 1000 )
endmacro (baci_test)

# Run test case for nested parallelism
macro (baci_test_Nested_Par arg1 arg2 restart)
  add_test(NAME ${arg1}-nestedPar
    COMMAND ${MPI_DIR}/bin/mpirun -np 3 $<TARGET_FILE:${baciname}> -ngroup=2 -glayout=1,2 -nptype=separateDatFiles ${PROJECT_SOURCE_DIR}/Input/${arg1}.dat xxx ${PROJECT_SOURCE_DIR}/Input/${arg2}.dat xxxAdditional)

  if (${restart})
    add_test(NAME ${arg1}-nestedPar-restart
    COMMAND ${MPI_DIR}/bin/mpirun -np 3 $<TARGET_FILE:${baciname}> -ngroup=2 -glayout=1,2 -nptype=separateDatFiles ${PROJECT_SOURCE_DIR}/Input/${arg1}.dat xxx restart=${restart} ${PROJECT_SOURCE_DIR}/Input/${arg2}.dat xxxAdditional restart=${restart})
  endif (${restart})

  set_tests_properties ( ${arg1}-nestedPar PROPERTIES TIMEOUT 200 )
endmacro (baci_test_Nested_Par)

# FRAMEWORK TESTS
macro (baci_framework_test arg nproc)
  set (RUNCUBIT ${CUBIT_DIR}/cubit\ -batch\ -nographics\ -nojournal\ ${PROJECT_SOURCE_DIR}/tests/framework-test/${arg}.jou)
  set (RUNPREEXODUS ./pre_exodus\ --exo=xxx_${arg}.e\ --bc=${PROJECT_SOURCE_DIR}/tests/framework-test/${arg}.bc\ --head=${PROJECT_SOURCE_DIR}/tests/framework-test/${arg}.head\ --dat=xxx.dat)
  set (RUNPRECURVE ./pre_curve\ --datfile=xxx.dat)
  set (RUNBACI ${MPI_DIR}/bin/mpirun\ -np\ ${nproc}\ $<TARGET_FILE:${baciname}>\ xxx.dat\ xxx)
  set (RUNPOSTFILTER ${MPI_DIR}/bin/mpirun\ -np\ ${nproc}\ ./post_drt_ensight\ --file=xxx)

  add_test(NAME ${arg}-p${nproc}-fw 
  COMMAND sh -c "${RUNCUBIT} && ${RUNPREEXODUS} && ${RUNPRECURVE} && ${RUNBACI} && ${RUNPOSTFILTER}")

# note: for the clean-up job in the end, every generated intermediate file has to start with "xxx"

  set_tests_properties ( ${arg}-p${nproc}-fw PROPERTIES TIMEOUT 1000 )
  set_tests_properties ( ${arg}-p${nproc}-fw PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR:; ERROR ;Error " )
  set_tests_properties ( ${arg}-p${nproc}-fw PROPERTIES ENVIRONMENT "PATH=$ENV{PATH}" )

endmacro (baci_framework_test)

# CUT TESTS
macro (cut_test nproc)
  set (RUNTESTS ${MPI_DIR}/bin/mpirun\ -np\ ${nproc}\ cut_test)
  add_test(NAME test-p${nproc}-cut 
  COMMAND sh -c "${RUNTESTS}")
  set_tests_properties ( test-p${nproc}-cut PROPERTIES TIMEOUT 1000 )
  set_tests_properties ( test-p${nproc}-cut PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR:; ERROR " )
endmacro (cut_test)

# POSTPROCESSING TEST
macro(post_processing arg nproc stresstype startstep)

  set (RUNPOSTFILTER_PAR ${MPI_DIR}/bin/mpirun\ -np\ ${nproc}\ ./post_drt_ensight\ --file=xxx\ --output=xxx_PAR\ --stress=${stresstype}\ --strain=${straintype}\ --start=${startstep})
  set (RUNPOSTFILTER_SER ./post_drt_ensight\ --file=xxx\ --output=xxx_SER\ --stress=${stresstype}\ --start=${startstep})
  add_test(NAME ${arg}-p${nproc}-pp  
  COMMAND sh -c " ${RUNPOSTFILTER_PAR} && ${RUNPOSTFILTER_SER} && pvpython\ ${PROJECT_SOURCE_DIR}/tests/post_processing_test/comparison.py")

  set_tests_properties ( ${arg}-p${nproc}-pp PROPERTIES TIMEOUT 1000 )
  set_tests_properties ( ${arg}-p${nproc}-pp PROPERTIES ENVIRONMENT "PATH=$ENV{PATH}" )
endmacro(post_processing)
 
###------------------------------------------------------------------ List of tests
##if(NOT TRILINOS_DEV)
baci_test(artaneurysm_hex 1 "")
baci_test(artaneurysm_hex 2 "")
post_processing(artaneurysm_hex 2 ndxyz 10)
baci_test(artaneurysm_tet 1 "")
baci_test(artaneurysm_tet 2 "")
post_processing(artaneurysm_tet 2 ndxyz 10)
baci_test(artaneurysm_hex_embed 1 "")
baci_test(artaneurysm_hex_embed 2 "")
baci_test(beam2rcubic 1 "")
baci_test(beam2rcubic 2 "")
baci_test(beam2rlinear 1 "")
baci_test(beam2rlinear 2 "")
baci_test(beam2rquadratic 1 "")
baci_test(beam2rquadratic 2 "")
baci_test(beam2rquartic 1 "")
baci_test(beam2rquartic 2 "")
baci_test(beam2test1 1 "")
baci_test(beam2test1 2 "")
baci_test(beam2test2 1 "")
baci_test(beam2test2 2 "")
baci_test(beam3contact_augmented 1 "")
baci_test(beam3contact_augmented 2 "")
baci_test(beam3contact_penalty 1 "")
baci_test(beam3contact_penalty 2 "")
baci_test(beam3contact_penalty_newgap_smoothed 1 "")
baci_test(beam3contact_penalty_newgap_smoothed 2 "")
baci_test(beam3eb_endmoment_quartercircle 1 "")
baci_test(beam3eb_endmoment_quartercircle 2 "")
baci_test(beam3ebtor_3Dendmoment 1 "")
baci_test(beam3ebtor_3Dendmoment 2 "")
baci_test(beam3ebtor_endmoment_endforce 1 "")
baci_test(beam3ebtor_endmoment_endforce 2 "")
baci_test(beam3iicontact_augmented 1 "")
baci_test(beam3iicontact_augmented 2 "")
baci_test(beam3iicontact_penalty 1 "")
baci_test(beam3iicontact_penalty 2 "")
baci_test(beam3contactstatmech_augmented 1 "")
baci_test(beam3contactstatmech_augmented 2 "")
baci_test(beam3cstatmechoctree_cobb 1 "")
baci_test(beam3cstatmechoctree_cobb 2 "")
baci_test(beam3cstatmechoctree_aabb 1 "")
baci_test(beam3cstatmechoctree_aabb 2 "")
baci_test(beam3contactstatmech_penalty 1 "")
baci_test(beam3contactstatmech_penalty 2 "")
baci_test(beam3iicontactstatmech_augmented 1 "")
baci_test(beam3iicontactstatmech_augmented 2 "")
baci_test(beam3iicontactstatmech_penalty 1 "")
baci_test(beam3iicontactstatmech_penalty 2 "")
baci_test(beam3test1cubic 1 "")
baci_test(beam3test1cubic 2 "")
baci_test(beam3test1linear 1 "")
baci_test(beam3test1linear 2 "")
baci_test(beam3test1quadratic 1 "")
baci_test(beam3test1quadratic 2 "")
baci_test(beam3test1quartic 1 "")
baci_test(beam3test1quartic 2 "")
baci_test(beam3test2cubic 1 "")
baci_test(beam3test2cubic 2 "")
baci_test(beam3test2linear 1 "")
baci_test(beam3test2linear 2 "")
baci_test(beam3test2quadratic 1 "")
baci_test(beam3test2quadratic 2 "")
baci_test(beam3test2quartic 1 "")
baci_test(beam3test2quartic 2 "")
baci_test(calcified_zugstab 2 "")
baci_test(combust_beltrami_8x8x8_ost 1 "")
baci_test(combust_beltrami_8x8x8_ost 2 "")
baci_test(combust_beltrami_8x8x8_ost_2iter 1 "")
baci_test(combust_beltrami_8x8x8_ost_2iter 2 "")
baci_test(combust_beltrami_8x8x8_afgenalpha_2iter 1 "")
baci_test(combust_beltrami_8x8x8_afgenalpha_2iter 2 "")
baci_test(combust_constpressurejump 1 "")
baci_test(combust_constpressurejump 2 "")
baci_test(combust_instat_timeInt_2Dduct 1 "")
baci_test(combust_instat_timeInt_2Dduct 2 "")
baci_test(combust_instat_timeInt_basic1Dduct 1 "")
baci_test(combust_instat_timeInt_basic1Dduct 2 3)
baci_test(combust_instat_timeIntEP_1Dduct_FGIFRS 1 "")
baci_test(combust_instat_timeIntEP_1Dduct_FGIFRS 2 "")
baci_test(combust_instat_timeIntEPNew_1Dduct 1 "")
baci_test(combust_instat_timeIntEPNew_1Dduct 2 "")
baci_test(combust_instat_timeIntMixed_1Dduct_FGIFRS_reverse 1 "")
baci_test(combust_instat_timeIntMixed_1Dduct_FGIFRS_reverse 2 "")
baci_test(combust_instat_timeIntSL_1Dduct_FGIFRS_reverse 1 "")
baci_test(combust_instat_timeIntSL_1Dduct_FGIFRS_reverse 2 "")
baci_test(combust_instat_timeIntSL_1Dduct_FGIFRS_redis 1 "")
baci_test(combust_instat_timeIntSL_1Dduct_FGIFRS_redis 2 "")
baci_test(combust_instat_timeIntSLEPNew_1Dduct 1 "")
baci_test(combust_instat_timeIntSLEPNew_1Dduct 2 "")
baci_test(combust_levelsetcut_horizontal 1 "")
baci_test(combust_levelsetcut_horizontal 2 "")
baci_test(combust_levelsetcut_vertical 1 "")
baci_test(combust_levelsetcut_vertical 2 "")
baci_test(combust_material_blending 1 "")
baci_test(combust_material_blending 2 "")
baci_test(combust_poiseuille 1 "")
baci_test(combust_poiseuille 2 "")
baci_test(combust_PDE_reinit_OST_instat_channel 1 "")
baci_test(combust_PDE_reinit_OST_instat_channel 2 "")
baci_test(combust_PDE_reinit_ICG_instat_channel 1 "")
baci_test(combust_PDE_reinit_ICG_instat_channel 2 "")
baci_test(combust_PDE_OST_circle_stat 1 "")
baci_test(combust_PDE_OST_circle_stat 2 "")
baci_test(combust_PDE_ICG_circle_stat 1 "")
baci_test(combust_PDE_ICG_circle_stat 2 "")
baci_test(condif1_straight_convection 1 "")
baci_test(condif1_straight_convection 2 "")
baci_test(condif2_skew_convection_10x10_maxent 1 "")
baci_test(condif2_straight_convection_direct 1 "")
baci_test(condif2_straight_convection_direct 2 "")
baci_test(condif2_straight_convection_vm3gls 2 "")
baci_test(condif3_straight_convection 1 "")
baci_test(condif3_straight_convection 2 "")
baci_test(constr2D_MPC_angle 1 "")
baci_test(constr2D_MPC_angle 2 "")
baci_test(constr2D_MPC_dist 1 "")
baci_test(constr2D_MPC_dist 2 "")
baci_test(constr2D_fsi 1 "")
baci_test(constr2D_fsi 2 "")
baci_test(constr2D_shrinking 1 "")
baci_test(constr2D_shrinking 2 "")
baci_test(constr3D_MPC_STI_trlike 1 "")
baci_test(constr3D_MPC_STI_trlike 2 "")
baci_test(constr3D_MPC_direct 1 "")
baci_test(constr3D_MPC_direct 2 "")
post_processing(constr3D_MPC_direct 2 ndxyz 25)
baci_test(constr3D_MPC_direct_sti 1 "")
baci_test(constr3D_MPC_direct_sti 2 "")
baci_test(constr3D_NormComp_STI 1 "")
baci_test(constr3D_NormComp_STI 2 "")
baci_test(constr3D_shrinking 1 "")
baci_test(constr3D_shrinking 2 "")
baci_test(constr3D_shrinking_ost_direct 1 "")
baci_test(constr3D_shrinking_ost_direct 2 "")
baci_test(constr3D_shrinking_ost_simple 1 "")
baci_test(constr3D_shrinking_ost_simple 2 "")
baci_test(constr3D_shrinking_ost_uzawa 1 "")
baci_test(constr3D_shrinking_ost_uzawa 2 "")
baci_test(constr3D_shrinking_stat_direct 1 "")
baci_test(constr3D_shrinking_stat_direct 2 "")
baci_test(constr3Dpen_shrinking_stat_direct 1 "")
baci_test(constr3Dpen_shrinking_stat_direct 2 "")
baci_test(contact2D_centrdiff 1 "")
baci_test(contact2D_centrdiff 2 "")
baci_test(contact2D_initfield 1 "")
baci_test(contact2D_initfield 2 "")
baci_test(contact2D_hertz_linstatic 1 "")
baci_test(contact2D_hertz_linstatic 2 "")
baci_test(contact2D_lin_augmented 1 "")
baci_test(contact2D_lin_augmented 2 "")
baci_test(contact2D_lin_duallagr 1 "")
baci_test(contact2D_lin_duallagr 2 "")
baci_test(contact2D_lin_penalty 1 "")
baci_test(contact2D_lin_penalty 2 "")
baci_test(contact2D_lin_petrovgalerkin 1 "")
baci_test(contact2D_lin_petrovgalerkin 2 "")
baci_test(contact2D_lin_stdlagr 1 "")
baci_test(contact2D_lin_stdlagr 2 "")
baci_test(contact2D_multibody 1 "")
baci_test(contact2D_multibody 2 "")
baci_test(contact2D_onlylin_stdlagr 1 "")
baci_test(contact2D_onlylin_stdlagr 2 "")
post_processing(contact2D_onlylin_stdlagr 2 ndxyz 50)
baci_test(contact2D_parredist 1 "")
baci_test(contact2D_parredist 2 "")
baci_test(contact2D_patch_linstatic 1 "")
baci_test(contact2D_patch_linstatic 2 "")
baci_test(contact2D_quad 1 "")
baci_test(contact2D_quad 2 "")
baci_test(contact2D_self 1 "")
baci_test(contact2D_self 2 "")
baci_test(contact2D_simpler 1 "")
baci_test(contact2D_simpler 2 "")
post_processing(contact2D_simpler 2 ndxyz 50)
baci_test(contact2D_slidingblock_coulomb 1 "")
baci_test(contact2D_slidingblock_coulomb 2 "")
baci_test(contact2D_slidingblock_tresca 1 "")
baci_test(contact2D_slidingblock_tresca 2 "")
baci_test(contact2D_twobeamshorizontal 1 "")
baci_test(contact2D_twobeamshorizontal 2 "")
baci_test(contact2D_twobeamsvertical 1 "")
baci_test(contact2D_twobeamsvertical 2 "")
baci_test(contact2D_wear 1 "")
baci_test(contact2D_wear 2 "")
baci_test(contact3D_lin_augmented 1 "")
baci_test(contact3D_lin_augmented 2 "")
baci_test(contact3D_lin_duallagr 1 "")
baci_test(contact3D_lin_duallagr 2 "")
baci_test(contact3D_lin_penalty 1 "")
baci_test(contact3D_lin_penalty 2 "")
baci_test(contact3D_lin_stdlagr 1 "")
baci_test(contact3D_lin_stdlagr 2 "")
baci_test(contact3D_onlylin_penalty_hex27 1 "")
baci_test(contact3D_onlylin_penalty_hex27 2 "")
baci_test(contact3D_patch_linstatic 1 "")
baci_test(contact3D_patch_linstatic 2 "")
baci_test(contact3D_pwlin_penalty_hex20 1 "")
baci_test(contact3D_pwlin_penalty_hex20 2 "")
post_processing(contact3D_pwlin_penalty_hex20 2 ndxyz 10)
baci_test(contact3D_pwlin_penalty_hex27 1 "")
baci_test(contact3D_pwlin_penalty_hex27 2 "")
post_processing(contact3D_pwlin_penalty_hex27 2 ndxyz 10)
baci_test(contact3D_pwlin_penalty_tet10 1 "")
baci_test(contact3D_pwlin_penalty_tet10 2 "")
baci_test(contact3D_quad_penalty_hex27 1 "")
baci_test(contact3D_quad_penalty_hex27 2 "")
baci_test(contact3D_quad_tet10 1 "")
baci_test(contact3D_quad_tet10 2 "")
baci_test(contact3D_simplebone1 2 "")
baci_test(contact3D_simplebone2 2 "")
baci_test(contact3D_slidingblock_augmented 1 "")
baci_test(contact3D_slidingblock_augmented 2 "")
baci_test(contact3D_slidingblock_augmented_tet10 1 "")
baci_test(contact3D_slidingblock_augmented_tet10 2 "")
baci_test(contact3D_slidingblock_lagrange_tet10 1 "")
baci_test(contact3D_slidingblock_lagrange_tet10 2 "")
baci_test(contact3D_slidingblock_lagrange_hex20 1 "")
baci_test(contact3D_slidingblock_lagrange_hex20 2 "")
baci_test(contact3D_slidingblock_coulomb 1 "")
baci_test(contact3D_slidingblock_coulomb 2 "")
baci_test(contact3D_slidingblock_penalty 1 "")
baci_test(contact3D_slidingblock_penalty 2 "")
baci_test(contact3D_wear 1 "")
baci_test(contact3D_wear 2 "")
baci_test(elch_1D_5ele_line3_3ions_stab_genalpha 1 "7")
baci_test(elch_1D_5ele_line3_3ions_stab_genalpha 2 "7")
baci_test(elch_1D_10ele_3ions_stab 1 "")
baci_test(elch_1D_10ele_3ions_stab 2 "")
baci_test(elch_1D_RDE_kineticsNewman 1 5)
baci_test(elch_1D_RDE_kineticsNewman 2 "")
baci_test(elch_1D_RDE_kineticsNewman_galvanostat 1 "10")
baci_test(elch_1D_RDE_kineticsNewman_2react 1 "")
baci_test(elch_1D_RDE_kineticsBard 1 "")
baci_test(elch_1D_RDE_kineticsBard 2 5)
baci_test(elch_1D_RDE_kineticsBard_galvanostat 1 "10")
baci_test(elch_1D_RDE_kineticsBard_2react 1 "")
baci_test(elch_1D_galvanostat 1 10)
baci_test(elch_Kwok_Wu 1 "")
baci_test(elch_Kwok_Wu 2 "")
baci_test(elch_Kwok_Wu_BDF2 2 "")
baci_test(elch_RCE_rotsymm_pbc 1 "")
baci_test(elch_RCE_rotsymm_pbc 2 "")
baci_test(elch_gaussian_hill_pbc_genalpha 2 8)
baci_test(elch_natconv 1 2)
baci_test(elch_natconv 2 2)
baci_test(elch_natconv_enc_pde_elim 1 2)
baci_test(elch_natconv_enc_pde_elim 2 2)
baci_test(elch_natconv_onescalar_boussinesq 1 "")
baci_test(elch_natconv_onescalar_boussinesq 2 "")
baci_test(elch_NonMatchingMesh_20x80 1 "")
baci_test(elch_NonMatchingMesh_20x80 2 "")
baci_test(elch_nurbs9 1 9)
baci_test(elch_nurbs9 2 9)
baci_test(elch_pper_butlervolmer 1 "")
baci_test(elch_pper_butlervolmer 2 "")
baci_test(elch_pper_primary_current_distribution 1 "")
baci_test(elch_pper_primary_current_distribution 2 "")
baci_test(elch_pper_secondary_current_distribution 1 "")
baci_test(elch_pper_secondary_current_distribution 2 "")
baci_test(elch_pulseplating_hex27 1 6)
baci_test(elch_pulseplating_hex27 2 6)
baci_test(elch_pulseplating_hex27_moving_boundary 2 "34")
baci_test(elch_straight_convection 1 "")
baci_test(elch_straight_convection 2 "")
baci_test(elch_test 1 "")
baci_test(elch_test 2 "")
baci_test(f2_channel20x20_drt_weak 1 "")
baci_test(f2_channel20x20_drt_weak 2 "")
baci_test(f2_channel_carreauyasuda_genalpha_drt 1 "")
baci_test(f2_channel_carreauyasuda_genalpha_drt 2 "")
baci_test(f2_channel_modpowerlaw_impl_drt 1 "")
baci_test(f2_channel_modpowerlaw_impl_drt 2 "")
baci_test(f2_dc20x20_drt_afg_cons 1 5)
baci_test(f2_dc20x20_drt_afg_cons 2 5)
baci_test(f2_dc20x20_drt_afgenalpha 1 5)
baci_test(f2_dc20x20_drt_afgenalpha 2 5)
baci_test(f2_drivencavity20x20_drt 1 5)
baci_test(f2_drivencavity20x20_drt 2 5)
baci_test(f2_drivencavity20x20_belos_drt 1 5)
baci_test(f2_drivencavity20x20_belos_drt 2 5)
baci_test(f2_drivencavity32x32_Re1000_stat 1 "")
baci_test(f2_drivencavity32x32_Re1000_stat 2 "")
baci_test(f2_drivencavity40x40_braess_sarazin_drt 1 "")
baci_test(f2_drivencavity40x40_braess_sarazin_drt 2 "")
baci_test(f2_freeosz20x20_drt 1 "")
baci_test(f2_freeosz20x20_drt 2 "")
baci_test(f2_impinging_flow_drt 1 "")
baci_test(f2_impinging_flow_drt 2 "")
baci_test(f2_nurbs9_dc_drt 1 "")
baci_test(f2_nurbs9_dc_drt 2 "")
baci_test(f2_ShearFlow_NonMatchingMesh_4x16_Con_Bmat 1 10)
baci_test(f2_ShearFlow_NonMatchingMesh_4x16_Con_Bmat 2 10)
baci_test(f2_ShearFlow_NonMatchingMesh_4x16_Con_Smat 1 "")
baci_test(f2_ShearFlow_NonMatchingMesh_4x16_Con_Smat 2 "")
baci_test(f2_ShearFlow_NonMatchingMesh_4x16_SP_Bmat 2 "")
baci_test(f2_ShearFlow_NonMatchingMesh_4x16_SP_Smat 1 "")
baci_test(f3_bel_tet4_drt_afg_cons 1 "")
baci_test(f3_bel_tet4_drt_afg_cons 2 "")
baci_test(f3_bel_tet4_drt_afgenalpha 1 "")
baci_test(f3_bel_tet4_drt_afgenalpha 2 "")
baci_test(f3_beltrami_tet4_drt 1 "")
baci_test(f3_beltrami_tet4_drt 2 "")
baci_test(f3_beltrami_8x8_afgenalpha_2iter_conv_cross_reynolds_gls 1 "")
baci_test(f3_beltrami_8x8_npgenalpha_2iter_conv_cross_reynolds_gls 1 "")
baci_test(f3_beltrami_8x8_genalpha_2iter_conv_cross_reynolds_gls 1 "")
baci_test(f3_beltrami_8x8_afgenalpha_066_2iter_conv_cross_reynolds_usfem 1 "")
baci_test(f3_beltrami_8x8_npgenalpha_066_2iter_conv_cross_reynolds_usfem 1 "")
baci_test(f3_beltrami_8x8_genalpha_066_2iter_conv_cross_reynolds_usfem 1 "")
baci_test(f3_beltrami_8x8_ost_066_2iter_conv_cross_reynolds_usfem 1 "")
baci_test(f3_beltrami_8x8x8_ost_conv 1 "")
baci_test(f3_beltrami_8x8x8_ost_2iter_conv 1 "")
baci_test(f3_beltrami_8x8x8_afgenalpha_2iter 1 "")
baci_test(f3_channel_carreauyasuda_genalpha_drt 1 "")
baci_test(f3_channel_carreauyasuda_genalpha_drt 2 "")
baci_test(f3_channel_modpowerlaw_impl_drt 1 "")
baci_test(f3_channel_modpowerlaw_impl_drt 2 "")
baci_test(f3_cylinder_001_ale_drt 1 "")
baci_test(f3_cylinder_001_ale_drt 2 "")
baci_test(f3_cylinder_001_drt 1 "")
baci_test(f3_cylinder_001_drt 2 "")
baci_test(f3_drivencavity_stat_drt 1 "")
baci_test(f3_drivencavity_stat_drt 2 "")
baci_test(f3_drivencavity_stat_belos_drt 1 "")
baci_test(f3_drivencavity_stat_belos_drt 2 "")
baci_test(f3_drt_genalpha_initialfield_volforce_pbc 1 "")
baci_test(f3_drt_genalpha_initialfield_volforce_pbc 2 "")
baci_test(f3_drt_impl_initialfield_volforce_pbc 1 "")
baci_test(f3_drt_impl_initialfield_volforce_pbc 2 "")
baci_test(f3_impl_channel_pressuregrad 1 "")
baci_test(f3_nurbs27_AfGenAlpha 3 "")
baci_test(f3_ew_Q1_27x_drt 1 "")
baci_test(f3_ew_Q1_27x_drt 2 "")
baci_test(f3_onewaycoupled_afgenalpha 2 8)
baci_test(f3_onewaycoupled_genalpha 2 8)
baci_test(f3_onewaycoupled_ost 2 8)
baci_test(f3_turbmodel_test_dynsmag 1 "")
baci_test(f3_turbmodel_test_dynsmag_dens 1 "")
baci_test(f3_turbmodel_test_fssgvisc 1 "")
baci_test(f3_turbmodel_test_mfs 1 "")
baci_test(f3_turbmodel_test_ssm 1 "")
baci_test(f3_turbmodel_test_ssmb 1 "")
baci_test(f3_turbinflow_doublechan 2 2)
baci_test(f3_turbinflow_doublechan_transferdbc 2 "")
baci_test(f3_simpler_dgs 1 "")
baci_test(f3_simpler_dgs 2 "")
baci_test(fs-kugel-tet 1 "")
baci_test(fs-kugel-tet 2 "")
baci_test(fs-quadrat 1 "")
baci_test(fs-quadrat 2 "")
baci_test(fs-wuerfel-tet 1 "")
baci_test(fs-wuerfel-tet 2 "")
baci_test(fs3i_part_1wc_finperm 1 10)
baci_test(fs3i_part_1wc_finperm 2 10)
baci_test(fs3i_part_1wc_infperm 1 10)
baci_test(fs3i_part_1wc_infperm 2 10)
baci_test(fsi_dc_mono_fs_ga_ga 1 8)
baci_test(fsi_dc_mono_fs_ga_ga 2 8)
baci_test(fsi_dc_mono_fs_ost_ga 1 "")
baci_test(fsi_dc_mono_fs_ost_ga 2 "")
baci_test(fsi_dc_mono_mtrfs_ost_ga 1 "")
baci_test(fsi_dc_mono_mtrfs_ost_ga 2 "")
baci_test(fsi_dc_mono_mtrss_ost_ga 1 "")
baci_test(fsi_dc_mono_mtrss_ost_ga 2 "")
baci_test(fsi_dc_mono_ss_ost_ga 1 "")
baci_test(fsi_dc_mono_ss_ost_ga 2 "")
baci_test(fsi_dc_part_ait_ga_ga 1 "")
baci_test(fsi_dc_part_ait_ga_ga 2 "")
baci_test(fsi_dc_part_ait_ost_ga 1 "")
baci_test(fsi_dc_part_ait_ost_ga 2 "")
baci_test(fsi_dc_part_mfnk_ost_ga 1 "")
baci_test(fsi_dc_part_mfnk_ost_ga 2 "")
baci_test(fsi_dc_part_mtrait_ost_ga 1 "")
baci_test(fsi_dc_part_mtrait_ost_ga 2 "")
baci_test(fsi_dc_part_sd_ost_ga 1 1)
baci_test(fsi_dc_part_sd_ost_ga 2 1)
baci_test(fsi_fp_mono_fs_bdf2_ga 1 8)
baci_test(fsi_fp_mono_fs_bdf2_ga 2 8)
baci_test(fsi_fp_mono_fs_bdf2_ost 1 8)
baci_test(fsi_fp_mono_fs_bdf2_ost 2 8)
baci_test(fsi_fp_mono_fs_ga_ga 1 8)
baci_test(fsi_fp_mono_fs_ga_ga 2 8)
baci_test(fsi_fp_mono_fs_ost_ga 1 8)
baci_test(fsi_fp_mono_fs_ost_ga 2 8)
baci_test(fsi_fp_mono_fs_ost_ost 1 8)
baci_test(fsi_fp_mono_fs_ost_ost 2 8)
baci_test(fsi_fp_mono_mtrfs_ga_ga 1 "")
baci_test(fsi_fp_mono_mtrfs_ga_ga 2 "")
baci_test(fsi_fp_mono_mtrfs_ost_ga 1 8)
baci_test(fsi_fp_mono_mtrfs_ost_ga 2 8)
baci_test(fsi_fp_mono_mtrss_ga_ga 1 8)
baci_test(fsi_fp_mono_mtrss_ga_ga 2 8)
baci_test(fsi_fp_mono_mtrss_ost_ga 1 8)
baci_test(fsi_fp_mono_mtrss_ost_ga 2 8)
baci_test(fsi_fp_mono_ss_bdf2_ga 1 8)
baci_test(fsi_fp_mono_ss_bdf2_ga 2 8)
baci_test(fsi_fp_mono_ss_bdf2_ost 1 8)
baci_test(fsi_fp_mono_ss_bdf2_ost 2 8)
baci_test(fsi_fp_mono_ss_ga_ga 1 8)
baci_test(fsi_fp_mono_ss_ga_ga 2 8)
baci_test(fsi_fp_mono_ss_ost_ga 1 8)
baci_test(fsi_fp_mono_ss_ost_ga 2 8)
baci_test(fsi_fp_mono_ss_ost_ost 1 8)
baci_test(fsi_fp_mono_ss_ost_ost 2 8)
baci_test(fsi_oval3D_sliding_up 1 "")
baci_test(fsi_oval3D_sliding_up 2 "")
baci_test(fsi_ow3D_mtr_drt 1 "")
baci_test(fsi_ow3D_mtr_drt 2 "")
baci_test(fsi_potential_LJ 1 "")
baci_test(fsi_potential_LJ_par 2 "")
baci_test(fsi_pw_mono_fs_ga_ga_contact 1 "")
baci_test(fsi_pw_mono_fs_ga_ga_contact 2 "")
baci_test(fsi_pw_mono_fs_ga_ga 1 3)
baci_test(fsi_pw_mono_fs_ga_ga 2 3)
baci_test(fsi_pw_mono_fs_ost_ga 1 3)
baci_test(fsi_pw_mono_fs_ost_ga 2 3)
baci_test(fsi_pw_mono_ss_ga_ga 1 "")
baci_test(fsi_pw_mono_ss_ga_ga 2 "")
baci_test(fsi_pw_mono_ss_ost_ga_fsiamg 1 3)
baci_test(fsi_pw_mono_ss_ost_ga_fsiamg 2 3)
baci_test(fsi_pw_mono_ss_ost_ga 1 3)
baci_test(fsi_pw_mono_ss_ost_ga 2 3)
baci_test(fsi_ring2D_sliding_tot 1 "")
baci_test(fsi_ring2D_sliding_tot 2 "")
baci_test(fsi_ring2D_sliding_up 1 "")
baci_test(fsi_ring2D_sliding_up 2 "")
baci_test(fsi_ves_mono_mtrfs_ost_ost_slidcurr 1 "")
baci_test(fsi_ves_mono_mtrfs_ost_ost_slidcurr 2 "")
baci_test(fsi_ves_mono_mtrfs_ost_ost_slidrot 1 "")
baci_test(fsi_ves_mono_mtrfs_ost_ost_slidrot 2 "")
baci_test(fsi_ves_mono_mtrss_ost_ost_slidcurr 1 "")
baci_test(fsi_ves_mono_mtrss_ost_ost_slidcurr 2 "")
baci_test(fsi_ves_mono_mtrss_ost_ost_slidrot 1 "")
baci_test(fsi_ves_mono_mtrss_ost_ost_slidrot 2 "")
baci_test(inv_analysis 1 "")
baci_test(invaaa 1 "")
baci_test(invaaa 2 "")
baci_test(inversedesign 1 "")
baci_test(inversedesign 2 "")
baci_test(inversedesign_tet4 1 "")
baci_test(inversedesign_tet4 2 "")
baci_test(loadlin 1 "")
baci_test(loadlin 2 "")
baci_test(locsys2D 1 "")
baci_test(locsys2D 2 "")
baci_test(locsys3D 1 "")
baci_test(locsys3D 2 "")
baci_test(locsys3D_sti 1 "")
baci_test(locsys3D_sti 2 "")
baci_test(loma_2d_heated_chan_30x100 1 10)
baci_test(loma_2d_heated_chan_30x100 2 10)
baci_test(loma_mono_2d_heated_chan_30x100 1 10)
baci_test(loma_mono_2d_heated_chan_30x100 2 10)
baci_test(loma_2d_heated_chan_sgs_30x100 1 10)
baci_test(loma_2d_heated_chan_sgs_30x100 2 10)
baci_test(loma_mono_2d_heated_chan_sgs_30x100 1 10)
baci_test(loma_mono_2d_heated_chan_sgs_30x100 2 10)
baci_test(loma_turbmodel_test_dynsmag 1 "")
baci_test(loma_turbmodel_test_fssgd 1 "")
baci_test(loma_turbmodel_test_mfs 1 "")
baci_test(mat_balzaniso_plate 1 "")
baci_test(mat_balzaniso_plate 2 "")
baci_test(mat_contchain_remo 1 "")
baci_test(mat_constraintmixture 1 "")
baci_test(mat_constraintmixture 2 "")
baci_test(mat_growth_ip 1 " ")
baci_test(mat_growth_ip 2 " ")
baci_test(mat_holzapfelcardio 1 "")
baci_test(mat_holzapfelcardio 2 "")
baci_test(mat_humphreycardio 1 "")
baci_test(mat_humphreycardio 2 "")
baci_test(meshtying2D_patch_linstatic 1 "")
baci_test(meshtying2D_patch_linstatic 2 "")
baci_test(meshtying2D_structure 1 "")
baci_test(meshtying2D_structure 2 "")
baci_test(meshtying2D_structure_crosspoints 1 "")
baci_test(meshtying2D_structure_crosspoints 2 "")
baci_test(meshtying3D_structure_augmented 1 "")
baci_test(meshtying3D_structure_augmented 2 "")
baci_test(meshtying3D_structure_duallagr 1 "")
baci_test(meshtying3D_structure_duallagr 2 "")
baci_test(meshtying3D_structure_errors 1 "")
baci_test(meshtying3D_structure_errors 2 "")
baci_test(meshtying3D_structure_penalty 1 "")
baci_test(meshtying3D_structure_penalty 2 "")
baci_test(meshtying3D_structure_simpler 1 "")
baci_test(meshtying3D_structure_simpler 2 "")
baci_test(meshtying3D_structure_stdlagr 1 "")
baci_test(meshtying3D_structure_stdlagr 2 "")
baci_test(meshtying3D_structure_tet10 1 "")
baci_test(meshtying3D_structure_tet10 2 "")
baci_test(mfsi-volconstr-fluidsplit 1 "")
baci_test(mfsi-volconstr-fluidsplit 2 "")
baci_test(mfsi-volconstr-structsplit 1 "")
baci_test(mfsi-volconstr-structsplit 2 "")
baci_test(mfsi_ow3D_mtr_drt 1 "")
baci_test(mfsi_ow3D_mtr_drt 2 "")
baci_test(mfsi_oval3D_sliding_up_fs 1 "")
baci_test(mfsi_oval3D_sliding_up_fs 2 "")
baci_test(orthopressure_shell8_eas_ans_sphere 1 8)
baci_test(orthopressure_shell8_eas_ans_sphere 2 8)
baci_test(orthopressure_solidh8_eas_sphere 1 "")
baci_test(orthopressure_solidh8_eas_sphere 2 "")
baci_test(orthopressure_solidh8_sphere 1 "")
baci_test(orthopressure_solidh8_sphere 2 "")
baci_test(orthopressure_solidw6_sphere 1 "")
baci_test(orthopressure_solidw6_sphere 2 "")
post_processing(orthopressure_solidw6_sphere 2 ndxyz 10)
baci_test(plastic_pressurisedcylinder 1 "")
baci_test(plastic_pressurisedcylinder 2 "")
baci_test(plastic_pressurisedcylinder_damage 1 "")
baci_test(plastic_pressurisedcylinder_damage 2 "")
baci_test(plastic_pressurisedcylinder_thrplastic 1 "")
baci_test(plastic_pressurisedcylinder_thrplastic 2 "")
baci_test(potential_surf_pointapprox_VDW_3D 1 "")
baci_test(potential_surf_pointapprox_VDW_3D 2 "")
baci_test(potential_vol_LJ_2D 1 "")
baci_test(potential_vol_LJ_2D 2 "")
baci_test(potential_vol_LJ_3D 1 "")
baci_test(potential_vol_LJ_3D 2 "")
baci_test(porous_cube5x5x5_compressed 1 "")
baci_test(porous_cube5x5x5_compressed 2 "")
baci_test(porous_cube5x5x5_constrained 1 "")
baci_test(porous_cube5x5x5_constrained 2 "")
baci_test(porous_cube5x5x5_comp_tet4 1 "")
baci_test(porous_cube5x5x5_comp_tet4 2 "")
baci_test(poro_fl_split_dbc 1 "")
baci_test(poro_fl_split_dbc 2 "")
baci_test(poro_fl_split_nbc 1 "")
baci_test(poro_fl_split_nbc 2 "")
baci_test(poro_scatra_normalflux 1 "")
baci_test(poro_scatra_normalflux 2 "")
baci_test(poro_st_split_dbc 1 "")
baci_test(poro_st_split_dbc 2 "")
baci_test(poro_st_split_nbc 1 "")
baci_test(poro_st_split_nbc 2 "")
baci_test(prestress_fsi 1 "")
baci_test(prestress_fsi 2 "")
baci_test(prestress_fsi_id 1 "")
baci_test(prestress_fsi_id 2 "")
baci_test(prestress_hex8 1 "")
baci_test(prestress_hex8 2 "")
baci_test(prestress_hex8fbar 1 "")
baci_test(prestress_hex8fbar 2 "")
baci_test(prestress_tet4 1 "")
baci_test(prestress_tet4 2 "")
baci_test(prestress_weg6 1 "")
baci_test(prestress_weg6 2 "")
baci_test(red_airway_tissue_coupling 1 "")
baci_test(red_airway_tissue_coupling 2 "")
baci_test(scatra_gaussian_hill_pbc 1 "")
baci_test(scatra_gaussian_hill_pbc 2 "")
baci_test(scatra_gaussian_hill_pbc_gena_cons 1 "")
baci_test(scatra_gaussian_hill_pbc_gena_cons 2 "")
baci_test(scatra_gaussian_hill_pbc_genalpha 1 "90")
baci_test(scatra_gaussian_hill_pbc_genalpha 2 "90")
baci_test(scatra_graetz_problem_system 1 5)
baci_test(scatra_graetz_problem_system 2 5)
baci_test(scatra_graetz_problem_system_ost 2 5)
baci_test(scatra_myocard_material 1 "")
baci_test(scatra_myocard_material 2 "")
baci_test(scatra_nurbs9_channel 1 "")
baci_test(scatra_nurbs9_channel 2 5)
baci_test(scatra_turbmodel_test_mfs 1 "")
baci_test(scatra_turbmodel_test_mfs_ml 1 "")
baci_test(scatra_turbmodel_test_dynsmag 1 "")
baci_test(so3_h27_patchtest_nl_cost_drt 1 "")
baci_test(so3_h27_patchtest_nl_cost_drt 2 "")
baci_test(sohex20_patchtest_nl_cost_drt 1 "")
baci_test(sohex20_patchtest_nl_cost_drt 2 "")
baci_test(sohex27_patchtest_nl_cost_drt 1 "")
baci_test(sohex27_patchtest_nl_cost_drt 2 "")
baci_test(sohex8_disp_patchtest 1 "")
baci_test(sohex8_disp_patchtest 2 "")
baci_test(sohex8_easfull_cooks_nl 1 "")
baci_test(sohex8_easfull_cooks_nl 2 "")
baci_test(sohex8_easfull_patchtest 1 "")
baci_test(sohex8_easfull_patchtest 2 "")
baci_test(sohex8_easmild_cooks_nl 1 "")
baci_test(sohex8_easmild_cooks_nl 2 "")
baci_test(sohex8_easmild_cooks_nl_dynsta 1 "")
baci_test(sohex8_easmild_cooks_nl_dynsta 2 "")
baci_test(sohex8_easmild_cooks_nl_dynstas 1 "")
baci_test(sohex8_easmild_cooks_nl_dynstas 2 "")
baci_test(sohex8_incompr_block 1 "")
baci_test(sohex8_incompr_block 2 "")
baci_test(sohex8_incompr_block_nln 1 "")
baci_test(sohex8_incompr_block_nln 2 "")
baci_test(sohex8_multiscale_macro 1 1)
baci_test(sohex8_multiscale_macro 2 1)
baci_test(sohex8_multiscale_macro_2micro 1 1)
baci_test(sohex8_multiscale_macro_2micro 2 1)
baci_test(sohex8_remodel 1 " ")
baci_test(sohex8_remodel 2 " ")
baci_test(sohex8_surfactant 1 30)
baci_test(sohex8_surfactant 2 30)
baci_test(sohex8fbar_cooks_nl 1 "")
baci_test(sohex8fbar_cooks_nl 2 "")
baci_test(sohex8p1j1_cooks_nln_drt 1 "")
baci_test(sohex8p1j1_cooks_nln_drt 2 "")
baci_test(solidh8_mat_neohooke_drt 1 "")
baci_test(solidh8_mat_neohooke_drt 2 "")
baci_test(sosh8_cooks_nl 1 "")
baci_test(sosh8_cooks_nl 2 "")
baci_test(sosh8_freeflying_ruler 1 "")
baci_test(sosh8_freeflying_ruler 2 "")
baci_test(sosh8_freeflying_ruler_sta 1 "")
baci_test(sosh8_freeflying_ruler_sta 2 "")
baci_test(sosh8_freeflying_ruler_sti 1 "")
baci_test(sosh8_freeflying_ruler_sti 2 "")
baci_test(sosh8_lin_bending_patchtest 1 "")
baci_test(sosh8_lin_bending_patchtest 2 "")
baci_test(sosh8_mooneyrivlin_rubbersheet 1 "")
baci_test(sosh8_mooneyrivlin_rubbersheet 2 "")
baci_test(sosh8_thin_nonlinear_pinched_hemisphere 1 "")
baci_test(sosh8_thin_nonlinear_pinched_hemisphere 2 "")
baci_test(sosh8p8_cooks_nl_incomp 1 "")
baci_test(sosh8p8_cooks_nl_incomp 2 "")
baci_test(sosh8p8_hybrid_patchtest 1 "")
baci_test(sosh8p8_hybrid_patchtest 2 "")
baci_test(sosh8p8_lin_bending_patchtest 1 "")
baci_test(sosh8p8_lin_bending_patchtest 2 "")
baci_test(sosh8p8_patchtest 1 "")
baci_test(sosh8p8_patchtest 2 "")
baci_test(soshw6_thin_nonlinear_pinched_hemisphere 1 "")
baci_test(soshw6_thin_nonlinear_pinched_hemisphere 2 "")
post_processing(soshw6_thin_nonlinear_pinched_hemisphere 2 ndxyz 20)
baci_test(sotet10_cooks_nl 1 "")
baci_test(sotet10_cooks_nl 2 "")
post_processing(sotet10_cooks_nl 2 ndxyz 100)
baci_test(sotet10_patchtest 1 "")
baci_test(sotet10_patchtest 2 "")
baci_test(sotet4_patchtest 1 "")
baci_test(sotet4_patchtest 2 "")
baci_test(soweg6_disp_patchtest 1 "")
baci_test(soweg6_disp_patchtest 2 "")
baci_test(stc-mltube-sosh8 1 "")
baci_test(stc-mltube-sosh8 2 "")
baci_test(stc-mltube-sosh8p8 1 "")
baci_test(stc-mltube-sosh8p8 2 "")
baci_test(stru_ale_wear 1 "")
baci_test(strufem_hex8_eas_cantilever 1 "")
baci_test(strufem_hex8_eas_cantilever 2 "")
baci_test(strufem_shell8_eas_cantilever 1 "")
baci_test(strufem_shell8_eas_cantilever 2 "")
baci_test(strufem_sosh8_eas_cantilever 1 "")
baci_test(strufem_sosh8_eas_cantilever 2 "")
baci_test(strufem_tet10_disp_cantilever 1 "")
baci_test(strufem_tet10_disp_cantilever 2 "")
post_processing(strufem_tet10_disp_cantilever 2 ndxyz 30)
baci_test(strufem_tet4_disp_cantilever 1 "")
baci_test(strufem_tet4_disp_cantilever 2 "")
post_processing(strufem_tet4_disp_cantilever 2 ndxyz 30)
baci_test(strufem_wedge6_disp_cantilever 1 "")
baci_test(strufem_wedge6_disp_cantilever 2 "")
post_processing(strufem_wedge6_disp_cantilever 2 ndxyz 30)
baci_test(thermo-line 1 "")
baci_test(thermo-line 2 "")
baci_test(thermo3D_FBC_ost 1 "")
baci_test(thermo3D_FBC_ost 2 "")
baci_test(thermo3D_FBC_statics 1 "")
baci_test(thermo3D_FBC_statics 2 "")
baci_test(thermoconvection_Tempn 1 "")
baci_test(thermoconvection_Tempn 2 "")
baci_test(thermoconvection_Tempnp 1 "")
baci_test(thermoconvection_Tempnp 2 "")
baci_test(thickcant_hex 1 "")
baci_test(thickcant_hex 2 "")
baci_test(thickcant_tet4 1 "")
baci_test(thickcant_tet4 2 "")
baci_test(topopt_adjoint_channel_pseudo_instat_const_vel_lin_pres 1 "")
baci_test(topopt_adjoint_channel_pseudo_instat_lin_vel_quad_pres 1 "")
baci_test(topopt_adjoint_channel_pseudo_instat_quad_vel_lin_pres 1 "")
baci_test(topopt_adjoint_channel_pseudo_instat_quad_vel_lin_pres_dens2_visc3 1 "")
baci_test(topopt_adjoint_channel_instat_primal_and_dual 1 "")
baci_test(topopt_adjoint_channel_instat_primal_and_dual 2 "")
baci_test(topopt_adjoint_channel_instat_quad_vel_lin_pres_dens2_visc3 1 "")
baci_test(topopt_adjoint_channel_instat_quad_vel_lin_pres_dens2_visc3 2 "")
baci_test(topopt_adjoint_channel_instat_varying_theta 1 "")
baci_test(topopt_adjoint_channel_instat_varying_theta 2 "")
baci_test(topopt_adjoint_channel_stat_const_vel_lin_pres 1 "")
baci_test(topopt_adjoint_channel_stat_lin_vel_quad_pres 1 "")
baci_test(topopt_adjoint_channel_stat_quad_vel_lin_pres 1 "")
baci_test(topopt_adjoint_channel_stat_quad_vel_lin_pres_dens2_visc3 1 "")
baci_test(topopt_adjoint_channel_stat_quad_vel_lin_pres_dens2_visc3 2 "")
baci_test(topopt_opti_snake_one_constr 1 "")
baci_test(topopt_opti_snake_one_constr 2 "")
baci_test(topopt_opti_snake_multiple_constr 1 "")
baci_test(topopt_opti_snake_multiple_constr 2 "")
baci_test(topopt_reaction_channel 1 "")
baci_test(topopt_reaction_channel 2 "")
baci_test(topopt_reaction_plate 1 "")
baci_test(topopt_reaction_plate 2 "")
baci_test(torsion2cosine 1 "")
baci_test(torsion2cosine 2 "")
baci_test(torsion2quadratic 1 "")
baci_test(torsion2quadratic 2 "")
baci_test(torsion3cosine 1 "")
baci_test(torsion3cosine 2 "")
baci_test(torsion3quadratic 1 "")
baci_test(torsion3quadratic 2 "")
baci_test(truss2test1 1 "")
baci_test(truss3test1 1 "")
baci_test(tsi_lindilatation_geolin 1 "")
baci_test(tsi_lincompression_1waydisp 1 "")
baci_test(tsi_lincompression_sequstaggdisp 1 "")
baci_test(tsi_lincompression_sequstaggtemp 1 "")
baci_test(tsi_lincompression_iterstaggdisp 1 "")
baci_test(tsi_lincompression_iterstaggtemp 1 "")
baci_test(tsi_lincompression_iterstaggaitken 1 "")
baci_test(tsi_lincompression_monolithic 1 "")
baci_test(tsi_heatflux_iterstagg 1 "")
baci_test(tsi_heatflux_iterstaggaitken 1 "")
baci_test(tsi_heatflux_monolithic 1 "")
baci_test(tsi_heatflux_flexoutsurf_monolithic 1 "")
baci_test(tsi_heatflux_flexoutsurf_monolithic 2 "")
baci_test(tsi_heatconvection_iterstagg 1 "")
baci_test(tsi_heatconvection_iterstagg 2 "")
baci_test(tsi_heatconvection_monolithic_fine 2 "")
baci_test(tsi_heatconvection_monolithic 1 "")
baci_test(tsi_heatconvection_monolithic 2 "")
baci_test(tsi_simplrocketnozzle_monolithic 1 "")
baci_test(tsi_simplrocketnozzle_monolithic 2 "")
baci_test(tsi_pressurisedcylinder_robinson 1 "")
baci_test(tsi_pressurisedcylinder_robinson 2 "")
baci_test(tsi_contact3D_twobodies 1 "")
baci_test(tsi_contact3D_twobodies 2 "")
baci_test(tsi_contact3D_heatproduction 1 "")
baci_test(tsi_contact3D_heatproduction 2 "")
baci_test(tsi_contact3D_twowaycoupling 1 "")
baci_test(tsi_contact3D_twowaycoupling 2 "")
baci_test(viscoAnisotropic_anisocreep 1 60)
baci_test(viscoAnisotropic_anisocreep 2 60)
baci_test(viscoAnisotropic_creep 1 60)
baci_test(viscoAnisotropic_creep 2 60)
baci_test(viscoNeoHooke_creep 1 60)
baci_test(viscoNeoHooke_creep 2 60)
baci_test(w1_dyn_dirich_drt 1 "")
baci_test(w1_dyn_dirich_drt 2 "")
baci_test(w1t6_patch_linear 1 "")
baci_test(w1t6_patch_linear 2 "")
post_processing(w1t6_patch_linear 2 ndxyz 1)
baci_test(w1q4_cantilever_eas 1 "")
baci_test(w1q4_cantilever_eas 2 "")
post_processing(w1q4_cantilever_eas 2 ndxyz 1)
baci_test(w1q4_patch_eas 1 "")
baci_test(w1q4_patch_eas 2 "")
post_processing(w1q4_patch_eas 2 ndxyz 1)
baci_test(w1q8_patch_nonlinear 1 "")
baci_test(w1q8_patch_nonlinear 2 "")
post_processing(w1q8_patch_nonlinear 2 ndxyz 1)
baci_test(w1q9_dyn_rot_drt_gemm 1 "")
baci_test(w1q9_dyn_rot_drt_gena 1 "")
baci_test(w1q9_dyn_rot_drt_gena 2 "")
post_processing(w1q9_dyn_rot_drt_gena 2 ndxyz 100)
baci_test(womersley_f3 2 49)
baci_test(womersley_f3_neumann_out 2 49)
baci_test(xfluid_channel_constWDBC_inclinedCut 1 "")
baci_test(xfluid_channel_constWDBC_inclinedCut 2 "")
baci_test(xfluid_channel_constWDBC_inclinedCutMomTes 2 "")
baci_test(xfluid_channel_constWDBC_inclinedCutMomMom 2 "")
baci_test(xfluid_channel_constWDBC_inclinedCutDiv 2 "")
baci_test(xfluid_channel_constWDBC_inclinedCut_bodyf_hex20 1 "")
baci_test(xfluid_channel_constWDBC_inclinedCut_bodyf_hex20 2 "")
baci_test(xfluid_channel_constWDBC_straightCut 1 "")
baci_test(xfluid_channel_constWDBC_straightCut 2 "")
baci_test(xfluid_channel_constWDBC_straightCutDiv 1 "")
baci_test(xfluid_channel_constWDBC_straightCutMomTes 1 "")
baci_test(xfluid_channel_constWDBC_straightCutMomMom 1 "")
baci_test(xfluid_channel_constWDBC_straightCut_hex20 1 "")
baci_test(xfluid_channel_constWDBC_straightCut_hex20 2 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut 1 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut 2 "")
baci_test(xfluid_channel_shearWDBC_inclinedCutMomTes 1 "")
baci_test(xfluid_channel_shearWDBC_inclinedCutMomMom 1 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyf 1 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyf 2 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyfMomTes 2 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyfMomMom 2 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyf_instat 1 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyf_instat 2 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyf_instatDiv 2 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyf_instatMomTes 1 "")
#baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyf_instatMomMom 1 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyf_Nitsche 1 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyf_NitscheMomTes 2 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyf_NitscheMomMom 2 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyf_NitscheDiv 1 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyf_hex20 1 "")
baci_test(xfluid_channel_shearWDBC_inclinedCut_bodyf_hex20 2 "")
baci_test(xfluid_channel_shearWDBC_inclineddoubleCut_bodyf 1 "")
baci_test(xfluid_channel_shearWDBC_inclineddoubleCut_bodyf 2 "")
baci_test(xfluid_fluid_pusher 1 "")
baci_test(xfluid_multibody_hex8 1 "")
baci_test(xfluid_multibody_hex8 2 "")
baci_test(xfluid_multibody_hex8_EOS_itsolver_MLprec_NIT 2 "")
baci_test(xfluid_thin_struct_hex8_EOS_itsolver_MLprec_NIT 2 "")
baci_test(xfluid_multibody_hex8_EOS_itsolver_MLprec_MSH 2 "")
baci_test(xfluid_stat_sliver_cut_instat 1 "")
baci_test(xfluid_stat_sliver_cut_instatMomTes 2 "")
baci_test(xfluid_stat_sliver_cut_instatMomMom 2 "")
baci_test(xfluid_stat_sliver_cut_instatDiv 2 "")
baci_test(xfluid_stat_sliver_cut_instat_EOS_GP 1 "")
baci_test(xfsi_fluid_pusher_instat_OST_1_Nitsche 1 "")
baci_test(xfsi_fluid_pusher_instat_OST_1 1 "")
baci_test(xfsi_fluid_pusher_instat_OST_1_Div 1 "")
baci_test(one_d_3_artery_network 1 9999)
baci_test(trachea 2 100)
baci_test(red_airway_acinus_of_3Exponential_alveolar_ducts 2 599)
baci_test(red_airway_2_acinus_of_3Exponential_alveolar_ducts 2 599)
baci_test(red_airway_2_acinus_of_Exponential_alveolar_duct_elems 2 599)
baci_test(red_airway_acinus_inter_acinar 2 890)
baci_test(beam3iitest1linear 1 "")
baci_test(beam3iitest1linear 2 "")
baci_test(beam3iitest2linear 1 "")
baci_test(beam3iitest2linear 2 "")
baci_test(inversedesign_wedge6 1 "")
baci_test(inversedesign_wedge6 2 "")
baci_test(xff_beltrami_TesTes 1 "")
#baci_test(xff_beltrami_MomTes 1 "")
baci_test(xff_beltrami_MomMom 1 "")
baci_test(xff_beltrami_Div 1 "")
baci_test(xff_channel_simple_stat 1 "")
#baci_test(xff_channel_simple_stat 2 "")
baci_test(xff_shear_flow_stat 1 "")
baci_test(xff_shear_flow_stat 2 "")
baci_test(xff_shear_flow_2D_stat_neumann 1 "")
baci_test(xff_shear_flow_2D_stat_neumann 2 "")
baci_test(xff_shear_flow_2D_stat_neumann_inflowstab_onesided 2 " ")
baci_test(xff_shear_flow_2D_stat_neumann_hex20 1 "") 
baci_test(xff_shear_flow_2D_stat_neumann_hex20 2 "") 
baci_test(xff_shear_flow_2D_instat_neumann 1 "")
baci_test(xff_shear_flow_2D_instat_neumann 2 "")
baci_test(xff_shear_flow_2D_instat_initialfield_OST 1 "")
baci_test(xff_shear_flow_2D_instat_initialfield_OST 2 "")
baci_test(xff_stat_cylinder_testcut 1 "")
baci_test(xff_stat_cylinder_testcut 3 "")
baci_test(xff_moving_circle_parainflow_OST 1 "")
baci_test(xff_moving_circle_parainflow_OST 3 "")
baci_test(xff_moving_circle_shearflow_OST 1 "")
baci_test(xff_moving_circle_shearflow_OST 3 "")
baci_test(xff_moving_circle_shearflow_OST_keepghostval 1 "")
baci_test(xff_moving_circle_shearflow_OST_keepghostval 3 "")
baci_test(xff_moving_circle_shearflow_OST_2 1 "")
baci_test(xff_moving_circle_shearflow_OST_2 3 "")
baci_test(xff_moving_circle_shearflow_moredofsets_OST 1 "")
baci_test(xff_moving_circle_shearflow_moredofsets_OST 3 "")
baci_test(xff_shear_flow_2D_stat_neumann_moredofset_quad 1 "")
baci_test(xff_shear_flow_2D_stat_neumann_moredofset_quad 2 "")
baci_test(xff_shear_flow_2D_stat_neum_moredofs_circ_Nitsche_alesided 1 "")
baci_test(xff_shear_flow_2D_stat_neum_moredofs_circ_Nitsche_alesided 2 "")
baci_test(xff_shear_flow_2D_stat_neum_moredofs_circ_Nitsche_alesided_inflowaveraged_velgradintstab 2 "")
baci_test(xff_shear_flow_2D_stat_neum_moredofs_circ_Nitsche_xsided 1 "")
baci_test(xff_shear_flow_2D_stat_neum_moredofs_circ_Nitsche_xsided 2 "")
baci_test(xff_shear_flow_2D_stat_neum_moredofs_circ_Nitsche_xsided_inflowstab_onesided 2 "")
baci_test(xff_shear_flow_2D_stat_neum_moredofs_circ_Nitsche_xsidedMomTes 2 "")
baci_test(xff_shear_flow_2D_stat_neum_moredofs_circ_Nitsche_xsidedMomMom 2 "")
baci_test(xff_shear_flow_2D_stat_neum_moredofs_circ_Nitsche_xsidedDiv 2 "")
baci_test(xff_shear_flow_2D_stat_neumann_moredofset_circle 1 "")
baci_test(xff_shear_flow_2D_stat_neumann_moredofset_circle 2 "")
baci_test(xffsi_monolithic_1D_channel_fixedALE 1 "")
baci_test(xffsi_monolithic_1D_channel_fsidirich_fixedAle 1 "")
baci_test(xffsi_monolithic_1D_channel_fullnewton 1 "")
baci_test(xffsi_monolithic_1D_channel_moreelements_fullnewton 2 "")
baci_test(xffsi_monolithic_stat_channel_fixedale_partitioned_edgebased 2 "")
baci_test(xffsi_monolithic_stat_channel_fixedale_partitioned_resbased 1 "")
baci_test(two_phase_beltrami_8x8x8_ost_2iter 1 "")
baci_test(two_phase_beltrami_8x8x8_ost_2iter 2 "")
baci_test(two_phase_beltrami_8x8x8_afgenalpha_2iter 1 "")
baci_test(two_phase_beltrami_8x8x8_afgenalpha_2iter 2 "")
baci_test(two_phase_bubble_variable_surftens 4 "")
baci_test(two_phase_couette_flow 1 "")
baci_test(two_phase_couette_flow 2 "")
baci_test(two_phase_couette_flow_vkink 1 "")
baci_test(two_phase_couette_flow_vkink 2 "")
baci_test(two_phase_couette_flow_none 1 "")
baci_test(two_phase_couette_flow_none 2 "")
baci_test(two_phase_rayleigh_taylor_instability_32x128 1 4)
baci_test(two_phase_rayleigh_taylor_instability_32x128 2 4)
baci_test(two_phase_rayleigh_taylor_instability_32x128_FGI 1 "")
baci_test(two_phase_rayleigh_taylor_instability_32x128_FGI 2 "")
baci_test(two_phase_rayleigh_taylor_instability_32x128_redis 3 "")
baci_test(two_phase_stat_bubble_surftens_20x20x20 2 "")
baci_test(two_phase_stat_bubble_surftens_20x20x20_pjump 2 "")
baci_test(two_phase_stat_bubble_surftens_50x50x1 1 "")
baci_test(two_phase_water_column_jumpsp_slp 1 "")
baci_test(two_phase_water_column_jumpsp_slp 2 "")
baci_test(two_phase_water_column_fem_sl 1 "")
baci_test(two_phase_water_column_fem_sl 2 "")
baci_test(two_phase_water_column_fem_don 1 "")
baci_test(two_phase_water_column_fem_don 2 "")
baci_test(two_phase_water_column_kinksp_donqs 1 "")
baci_test(two_phase_water_column_kinksp_donqs 2 "")
baci_test(art_aorta_2D_lognorm_field 1 "")
baci_test(art_aorta_3D_lognorm_field_FFT 1 "")
baci_test_Nested_Par(tsi_heatconvection_monolithic tsi_heatconvection_monolithic "")
baci_test_Nested_Par(sohex8_multiscale_macro sohex8_multiscale_npsupport "1")
baci_test_Nested_Par(sohex8_multiscale_macro_2micro sohex8_multiscale_npsupport "1")

# testing the whole framework: from cubit via pre_exodus and baci to the post-filter:
baci_framework_test(tutorial_fluid 2)
baci_framework_test(tutorial_fsi 2)

# cut test 
cut_test(2)

if (HAVE_MueLu)
baci_test(f2_drivencavity20x20_muelu_drt 1 "")
baci_test(f2_drivencavity20x20_muelu_drt 2 "")
endif (HAVE_MueLu)

if (FS3I_INCACouplingTest)
  include (/home/hammerl/buildslave/FS3I_INCA_coupling_test.cmake)
endif (FS3I_INCACouplingTest)

# remove any output files from our tests
add_test(NAME test_cleanup COMMAND sh -c "rm -vf xxx* core.")





