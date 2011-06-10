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
      add_definitions("-DCHANGEDDATE=${baci_WC_LAST_CHANGED_DATE}")
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
add_definitions("-DPARALLEL")
add_definitions("-DPARMETIS")
add_definitions("-DLOOPL=15")
add_definitions("-DCCADISCRET")
add_definitions("-DTRILINOS_PACKAGE")
add_definitions("-DD_FLUID")
add_definitions("-DD_FSI")

# external library variable
set(LIBRARIES "")

#
### compile time options

option (DEBUG "Turn on special debug code" OFF)
if (DEBUG)
  add_definitions("-DDEBUG")
endif (DEBUG)

option (DEBUGCUTLIBRARY "Turn on special output inside the cut library to be able to track errors" OFF)
if (DEBUGCUTLIBRARY)
  add_definitions("-DDEBUGCUTLIBRARY")
endif (DEBUGCUTLIBRARY)

option (STKADAPTIVE "adaptive meshes based on the stk library" OFF)
if (STKADAPTIVE)
  add_definitions("-DSTKADAPTIVE")
endif (STKADAPTIVE)

#   additional features

option (BINIO "binary io" ON)
if (BINIO)
  find_package(HDF5 REQUIRED)
  if (HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
    link_directories(${HDF5_LIBRARY_DIRS})
    set(LIBRARIES
      ${LIBRARIES}
      ${HDF5_LIBRARIES}
      )
  endif(HDF5_FOUND)
  add_definitions("-DBINIO")
endif(BINIO)

find_package(Qhull REQUIRED)
if (QHULL_FOUND)
  include_directories(${QHULL_INCLUDE_DIRS})
  link_directories(${QHULL_LIBRARY_DIRS})
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

if (NOT INCLUDE_INSTALL_DIR)
  set ( INCLUDE_INSTALL_DIR "/lnm/lib/fc8/openmpi/include" )
endif (NOT INCLUDE_INSTALL_DIR)

if (NOT LIB_INSTALL_DIR)
  set ( LIB_INSTALL_DIR "/lnm/lib/fc8/openmpi/lib" )
endif (NOT LIB_INSTALL_DIR)

if (NOT Trilinos_DIR)
  set ( Trilinos_DIR "/lnm/lib/fc8/openmpi/trilinos" )
endif (NOT Trilinos_DIR)

#
### External Libraries

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
#find_package(MPI REQUIRED)
find_package(Trilinos REQUIRED)

#find_package(Blitz)

#set(Boost_USE_STATIC_LIBS ON) 
#set(Boost_DETAILED_FAILURE_MSG ON)
#set(Boost_USE_MULTITHREAD OFF) 
#set(Boost_DEBUG ON)
#FIND_PACKAGE( Boost COMPONENTS date_time filesystem iostreams program_options )
find_package( Boost COMPONENTS graph program_options )

if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  set(LIBRARIES
    ${Boost_LIBRARIES}
    ${LIBRARIES}
    )
else()
 message(FATAL_ERROR ${Boost_ERROR_REASON})
endif()

# trilinos libraries
# Choose trilinos libraries explicitly.
if (Trilinos_FOUND)
  include_directories(${Trilinos_INCLUDE_DIRS})
  include_directories(${Trilinos_TPL_INCLUDE_DIRS})
  link_directories(${Trilinos_LIBRARY_DIRS})
#  set(USED_Trilinos_LIBRARIES
#    trilinoscouplings
#    noxepetra
#    nox
#    ifpack
#    aztecoo
#    amesos
#    ml
#    zoltan
#    epetraext
#    epetra
#    teuchos
#    )
  
#  # some parts of baci might require additional trilinos libraries

#  if (STKADAPTIVE)
#    set(USED_Trilinos_LIBRARIES
#      stk_search stk_linsys
#      stk_mesh_fem stk_mesh_base 
#      #stk_util_use_cases 
#      stk_util_parallel 
#      #stk_util_diag 
#      #stk_util_env 
#      stk_util_util
#      #intrepid
#      #fei_test_utils fei_trilinos fei_base
#      shards
#      #isorropia
#      #zoltan
#      ${USED_Trilinos_LIBRARIES}
#      )
#  endif (STKADAPTIVE)

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

#find_package(Umfpack)
#if (UMFPACK_FOUND)
#  include_directories(${UMFPACK_INCLUDE_DIRS})
#  link_directories(${UMFPACK_LIBRARY_DIRS})
#  set(LIBRARIES
#    ${LIBRARIES}
#    ${UMFPACK_LIBRARIES}
#    )
#else (UMFPACK_FOUND)
#  add_definitions("-DHAVENOT_UMFPACK")
#endif (UMFPACK_FOUND)

#find_package(SuperLU)
#if (SUPERLU_FOUND)
#  include_directories(${SUPERLU_INCLUDE_DIRS})
#  link_directories(${SUPERLU_LIBRARY_DIRS})
#  set(LIBRARIES
#    ${LIBRARIES}
#    ${SUPERLU_LIBRARIES}
#    )
#else (SUPERLU_FOUND)
#  add_definitions("-DHAVENOT_SUPERLU")
#endif (SUPERLU_FOUND)

#find_package(PARMETIS REQUIRED)
#if (PARMETIS_FOUND)
#  include_directories(${PARMETIS_INCLUDE_DIRS})
#  link_directories(${PARMETIS_LIBRARY_DIRS})
#  set(LIBRARIES
#    ${LIBRARIES}
#    ${PARMETIS_LIBRARIES}
#    )
#endif (PARMETIS_FOUND)

#if (Boost_FOUND)
#  include_directories(${Boost_INCLUDE_DIRS})
#  link_directories(${Boost_LIBRARY_DIRS})
#  set(LIBRARIES
#    ${LIBRARIES}
#    ${Boost_LIBRARIES}
#    )
#endif(Boost_FOUND)

#find_package(NemesisI REQUIRED)
#if (NEMESIS_FOUND)
#  include_directories(${NEMESIS_INCLUDE_DIRS})
#  link_directories(${NEMESIS_LIBRARY_DIRS})
#  set(LIBRARIES
#    ${LIBRARIES}
#    ${NEMESIS_LIBRARIES}
#    )
#endif(NEMESIS_FOUND)

#find_package(ExodusII REQUIRED)
#if (EXODUS_FOUND)
#  include_directories(${EXODUS_INCLUDE_DIRS})
#  link_directories(${EXODUS_LIBRARY_DIRS})
#  set(LIBRARIES
#    ${LIBRARIES}
#    ${EXODUS_LIBRARIES}
#    )
#endif(EXODUS_FOUND)

find_package(Netcdf REQUIRED)
if (NETCDF_FOUND)
  include_directories(${NETCDF_INCLUDE_DIRS})
  link_directories(${NETCDF_LIBRARY_DIRS})
  set(LIBRARIES
    ${LIBRARIES}
    ${NETCDF_LIBRARIES}
    )
endif(NETCDF_FOUND)

#if (LAPACK_FOUND)
#  include_directories(${LAPACK_INCLUDE_DIRS})
#  link_directories(${LAPACK_LIBRARY_DIRS})
#  set(LIBRARIES
#    ${LIBRARIES}
#    ${LAPACK_LIBRARIES}
#    )
#endif(LAPACK_FOUND)

#if (BLAS_FOUND)
#  include_directories(${BLAS_INCLUDE_DIRS})
#  link_directories(${BLAS_LIBRARY_DIRS})
#  set(LIBRARIES
#    ${LIBRARIES}
#    ${BLAS_LIBRARIES}
#    )
#endif(BLAS_FOUND)

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
  src/drt_fluid/drt_pbcdofset.cpp
  src/drt_fluid/drt_periodicbc.cpp
  src/drt_lib/drt_assemblestrategy.cpp
  src/drt_lib/drt_condition.cpp
  src/drt_lib/drt_conditiondefinition.cpp
  src/drt_lib/drt_container.cpp
  src/drt_lib/drt_discret.cpp
  src/drt_lib/drt_discret_conditions.cpp
  src/drt_lib/drt_discret_evaluate.cpp
  src/drt_lib/drt_discret_fillcomplete.cpp
  src/drt_lib/drt_discret_partition.cpp
  src/drt_lib/drt_discret_utils.cpp
  src/drt_lib/drt_dofset.cpp
  src/drt_lib/drt_dofset_base.cpp
  src/drt_lib/drt_dofset_proxy.cpp
  src/drt_lib/drt_dofset_transparent.cpp
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
  src/drt_lib/drt_utils_nullspace.cpp
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
  src/drt_mat/carreauyasuda.cpp
  src/drt_mat/charmm.cpp
  src/drt_mat/cnst_1d_art.cpp
  src/drt_mat/compogden.cpp
  src/drt_mat/constraintmixture.cpp
  src/drt_mat/contchainnetw.cpp
  src/drt_mat/elasthyper.cpp
  src/drt_mat/ferech_pv.cpp
  src/drt_mat/fourieriso.cpp
  src/drt_mat/growth_ip.cpp
  src/drt_mat/holzapfelcardiovascular.cpp
  src/drt_mat/humphreycardiovascular.cpp
  src/drt_mat/ion.cpp
  src/drt_mat/itskov.cpp
  src/drt_mat/logneohooke.cpp
  src/drt_mat/lung_ogden.cpp
  src/drt_mat/lung_penalty.cpp
  src/drt_mat/material.cpp
  src/drt_mat/matlist.cpp
  src/drt_mat/matpar_bundle.cpp
  src/drt_mat/matpar_material.cpp
  src/drt_mat/matpar_parameter.cpp
  src/drt_mat/micromaterial.cpp
  src/drt_mat/mixfrac.cpp
  src/drt_mat/modpowerlaw.cpp
  src/drt_mat/mooneyrivlin.cpp
  src/drt_mat/neohooke.cpp
  src/drt_mat/newtonianfluid.cpp
  src/drt_mat/permeablefluid.cpp
  src/drt_mat/plasticneohooke.cpp
  src/drt_mat/plasticlinelast.cpp
  src/drt_mat/protein.cpp
  src/drt_mat/scatra_mat.cpp
  src/drt_mat/stvenantkirchhoff.cpp
  src/drt_mat/sutherland.cpp
  src/drt_mat/thermostvenantkirchhoff.cpp
  src/drt_mat/viscoanisotropic.cpp
  src/drt_mat/visconeohooke.cpp
  src/drt_mat/yeoh.cpp
  src/drt_matelast/elast_coupanisoexpotwo.cpp
  src/drt_matelast/elast_coupanisoneohooketwo.cpp
  src/drt_matelast/elast_coupblatzko.cpp
  src/drt_matelast/elast_couplogneohooke.cpp
  src/drt_matelast/elast_coupvarga.cpp
  src/drt_matelast/elast_isocub.cpp
  src/drt_matelast/elast_isoexpo.cpp
  src/drt_matelast/elast_isomooneyrivlin.cpp
  src/drt_matelast/elast_isoneohooke.cpp
  src/drt_matelast/elast_isoquad.cpp
  src/drt_matelast/elast_isovarga.cpp
  src/drt_matelast/elast_isoyeoh.cpp
  src/drt_matelast/elast_summand.cpp
  src/drt_matelast/elast_varisoneohooke.cpp
  src/drt_matelast/elast_varisoquad.cpp
  src/drt_matelast/elast_varisocub.cpp
  src/drt_matelast/elast_varisoexpo.cpp
  src/drt_matelast/elast_vologden.cpp
  src/drt_matelast/elast_volpenalty.cpp
  src/drt_matelast/elast_volsussmanbathe.cpp
  src/drt_patspec/patspec.cpp
)
set(OBJS_DRT_STATMECH
  src/drt_statmech/brownianrod.cpp
  src/drt_statmech/statmech_manager.cpp
  src/drt_statmech/statmech_output.cpp
  src/drt_statmech/statmech_time.cpp
)
set(OBJS_DRT_ADAPT
  src/drt_adapter/adapter_algorithmbase.cpp
  src/drt_adapter/adapter_coupling.cpp
  src/drt_adapter/adapter_coupling_mortar.cpp
  src/drt_adapter/adapter_fluid.cpp
  src/drt_adapter/adapter_fluid_ale.cpp
  src/drt_adapter/adapter_fluid_base_algorithm.cpp
  src/drt_adapter/adapter_fluid_combust.cpp
  src/drt_adapter/adapter_fluid_genalpha.cpp
  src/drt_adapter/adapter_fluid_impl.cpp
  src/drt_adapter/adapter_fluid_lung.cpp
  src/drt_adapter/adapter_fluid_moving_boundary.cpp
  src/drt_adapter/adapter_fluid_projection.cpp
  src/drt_adapter/adapter_fluid_xfem.cpp
  src/drt_adapter/adapter_scatra_base_algorithm.cpp
  src/drt_adapter/adapter_scatra_fluid_ale_coupling_algo.cpp
  src/drt_adapter/adapter_scatra_fluid_coupling_algorithm.cpp
  src/drt_adapter/adapter_structure.cpp
  src/drt_adapter/adapter_structure_constr_merged.cpp
  src/drt_adapter/adapter_structure_lung.cpp
  src/drt_adapter/adapter_structure_strugenalpha.cpp
  src/drt_adapter/adapter_structure_timint.cpp
  src/drt_adapter/adapter_structure_timint_adaptive.cpp
  src/drt_adapter/adapter_structure_timint_expl.cpp
  src/drt_adapter/adapter_structure_wrapper.cpp
  src/drt_adapter/adapter_thermo.cpp
  src/drt_adapter/adapter_thermo_timint.cpp
  src/drt_adapter/adapter_xfluid_impl.cpp
  src/drt_adapter/adapter_fluid_xfluid_impl.cpp	
  src/drt_adapter/adapter_fluid_fluid_impl.cpp
)
set(OBJS_DRT_CONTACT_LIB
  src/drt_contact/contact_element.cpp
  src/drt_contact/contact_node.cpp
  src/drt_contact/friction_node.cpp
)
set(OBJS_DRT_RED_AIRWAYS_LIB
  src/drt_red_airways/airway_evaluate.cpp
  src/drt_red_airways/airway_impl.cpp
  src/drt_red_airways/red_airway.cpp
  src/drt_red_airways/red_airway_input.cpp
  src/drt_red_airways/red_airway_resulttest.cpp
)
set(OBJS_DRT_INPAR
  src/drt_inpar/drt_boolifyparameters.cpp
  src/drt_inpar/drt_validconditions.cpp
  src/drt_inpar/drt_validmaterials.cpp
  src/drt_inpar/drt_validparameters.cpp
)
set(OBJS_DRT_ARTNET
  src/drt_art_net/art_net_dyn_drt.cpp
  src/drt_art_net/artnetexplicitintegration.cpp
)
set(OBJS_INTERSECTION
  src/drt_xfem/dof_distribution_switcher.cpp
  src/drt_xfem/dof_management.cpp
  src/drt_xfem/dof_management_element.cpp
  src/drt_xfem/enrichment.cpp
  src/drt_xfem/enrichmentvalues.cpp
  src/drt_xfem/enrichment_utils.cpp
  src/drt_xfem/field_enriched.cpp
  src/drt_xfem/interface.cpp
  src/drt_xfem/interfacexfsi.cpp
  src/drt_xfem/physics.cpp
  src/drt_xfem/spacetime_boundary.cpp
  src/drt_xfem/startvalues.cpp
  src/drt_xfem/xdofmapcreation.cpp
  src/drt_xfem/xdofmapcreation_combust.cpp
  src/drt_xfem/xdofmapcreation_fsi.cpp
  src/drt_xfem/xdofmapcreation_parallel_utils.cpp
  src/drt_xfem/xfem_fluiddofset.cpp
  src/drt_xfem/xfem_fluidwizard.cpp
  src/drt_geometry/geo_intersection.cpp
  src/drt_combust/combust_interface.cpp
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
)
set(OBJS_POST_DRT_Q9NECONVERT
  src/post_q9neconvert/post_drt_q9neconvert.cpp
)
set(OBJS_DRT_COMBUST
  src/drt_combust/combust_algorithm.cpp
  src/drt_combust/combust_dyn.cpp
  src/drt_combust/combust_flamefront.cpp
  src/drt_combust/combust_fluidimplicitintegration.cpp
  src/drt_combust/combust_fluidresulttest.cpp
#  src/drt_combust/combust_interface.cpp
  src/drt_combust/combust_refinementcell.cpp
  src/drt_combust/combust_reinitializer.cpp
  src/drt_combust/combust_utils.cpp
)
set(OBJS_DRT_BEAMCONTACT
  src/drt_beamcontact/beam3contactstrugenalpha.cpp
)
set(OBJS_FILTER_COMMON
  src/filter_common/filter_evaluation.cpp
)
set(OBJS_DRT_BEAMCONTACT_LIB
  src/drt_beamcontact/beam3contact.cpp
  src/drt_beamcontact/beam3contact_manager.cpp
  src/drt_beamcontact/beamcontact_utilitis.cpp
  src/drt_beamcontact/beam3contact_octtree.cpp
)
set(OBJS_DRT_STRU_MULTI
  src/drt_stru_multi/convert_mat.cpp
  src/drt_stru_multi/microstatic.cpp
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
  src/drt_so3/so_hex8_thermo_lin_evaluate.cpp
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
  src/drt_so3/so_nstet.cpp
  src/drt_so3/so_nstet_evaluate.cpp
  src/drt_so3/so_nstet_input.cpp
  src/drt_so3/so_nstet_nodalstrain.cpp
  src/drt_so3/so_nurbs27.cpp
  src/drt_so3/so_nurbs27_evaluate.cpp
  src/drt_so3/so_nurbs27_input.cpp
  src/drt_so3/so_ptet.cpp
  src/drt_so3/so_ptet_evaluate.cpp
  src/drt_so3/so_ptet_input.cpp
  src/drt_so3/so_ptet_nodalstrain.cpp
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
  src/drt_so3/so_tet4_input.cpp
  src/drt_so3/so_weg6.cpp
  src/drt_so3/so_weg6_evaluate.cpp
  src/drt_so3/so_weg6_input.cpp
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
set(OBJS_DRT_UTILS
  src/drt_fem_general/debug_nan.cpp
  src/drt_fem_general/drt_utils_boundary_integration.cpp
  src/drt_fem_general/drt_utils_dgfem_basisfunctions.cpp
  src/drt_fem_general/drt_utils_fem_shapefunctions.cpp
  src/drt_fem_general/drt_utils_integration.cpp
  src/drt_fem_general/drt_utils_gausspoints.cpp
  src/drt_fem_general/drt_utils_local_connectivity_matrices.cpp
  src/drt_fem_general/largerotations.cpp
  src/drt_fluid/fluid_rotsym_periodicbc_utils.cpp
  src/drt_geometry/integrationcell.cpp
  src/drt_geometry/intersection.cpp
  src/drt_geometry/intersection_interfacepoint.cpp
  src/drt_geometry/intersection_math.cpp
  src/drt_geometry/intersection_service.cpp
  src/drt_geometry/position_array.cpp
  src/drt_geometry/searchtree.cpp
  src/drt_geometry/searchtree_geometry_service.cpp
  src/drt_geometry/searchtree_nearestobject.cpp
  src/drt_geometry/tetrahedradecomposition.cpp
  src/drt_io/io_gmsh.cpp
  src/drt_io/io_gmsh_xfem_extension.cpp
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
)
set(OBJS_DRT_TSI
  src/drt_tsi/tsi_algorithm.cpp
  src/drt_tsi/tsi_monolithic.cpp
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
  src/math/math1.c
  src/shell8/s8_mat_linel.c
  src/shell8/s8_mat_ogden.c
  src/shell8/s8_mat_ogden2.c
  src/shell8/s8_mat_ogden3.c
  src/shell8/s8_mattransform.c
  src/shell8/s8_tmat.c
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
  src/drt_fsi/fsi_lung_scatra.cpp
  src/drt_fsi/fsi_fluid_ale.cpp
  src/drt_fsi/fsi_fluid_xfem.cpp
  src/drt_fsi/fsi_lagrangeprec.cpp
  src/drt_fsi/fsi_lungmonolithic.cpp
  src/drt_fsi/fsi_lungmonolithic_fluidsplit.cpp
  src/drt_fsi/fsi_lungmonolithic_structuresplit.cpp
  src/drt_fsi/fsi_matrixtransform.cpp
  src/drt_fsi/fsi_monolithic.cpp
  src/drt_fsi/fsi_monolithic_linearsystem.cpp
  src/drt_fsi/fsi_monolithiclagrange.cpp
  src/drt_fsi/fsi_monolithicfluidsplit.cpp
  src/drt_fsi/fsi_monolithicstructuresplit.cpp
  src/drt_fsi/fsi_monolithicxfem.cpp
  src/drt_fsi/fsi_mortarmonolithic_structuresplit.cpp
  src/drt_fsi/fsi_mortarmonolithic_fluidsplit.cpp
  src/drt_fsi/fsi_nox_aitken.cpp
  src/drt_fsi/fsi_nox_epsilon.cpp
  src/drt_fsi/fsi_nox_fixpoint.cpp
  src/drt_fsi/fsi_nox_group.cpp
  src/drt_fsi/fsi_nox_jacobian.cpp
  src/drt_fsi/fsi_nox_linearsystem_bgs.cpp
  src/drt_fsi/fsi_nox_linearsystem_gcr.cpp
  src/drt_fsi/fsi_nox_linearsystem_partitioned.cpp
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
  src/drt_fsi/fsi_partitionedmonolithic.cpp
  src/drt_fsi/fsi_robin.cpp
  src/drt_fsi/fsi_robinneumann.cpp
  src/drt_fsi/fsi_statustest.cpp
  src/drt_fsi/fsi_structureale.cpp
  src/drt_fsi/fsi_utils.cpp
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
  src/drt_scatra/scatra_ele_impl_utils.cpp
  src/drt_scatra/scatra_element.cpp
  src/drt_scatra/scatra_element_boundary_evaluate.cpp
  src/drt_scatra/scatra_element_evaluate.cpp
  src/drt_scatra/scatra_element_input.cpp
)
set(OBJS_DRT_T2_LIB
  src/drt_truss2/truss2.cpp
  src/drt_truss2/truss2_evaluate.cpp
  src/drt_truss2/truss2_input.cpp
)
set(OBJS_DRT_STATMECH_LIB
  src/drt_statmech/bromotion_manager.cpp
)
set(OBJS_DRT_COMBUST3_LIB
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
  src/drt_thermo/thrtimint_impl_contact.cpp
  src/drt_thermo/thrtimint_ost.cpp
  src/drt_thermo/thrtimint_statics.cpp
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
  src/drt_fluid/fluid_projectionmethod.cpp
  src/drt_fluid/fluid_meshtying.cpp
  src/drt_fluid/fluid_utils.cpp
  src/drt_fluid/fluid_utils_mapextractor.cpp
  src/drt_fluid/fluid_MHD_evaluate.cpp
  src/drt_fluid/fluid_windkessel_optimization.cpp
  src/drt_fluid/fluidimpedancecondition.cpp
  src/drt_fluid/fluid_volumetric_surfaceFlow_condition.cpp
  src/drt_fluid/fluidimplicitintegration.cpp
  src/drt_fluid/fluidxfluidimplicitintegration.cpp
  src/drt_fluid/fluidxfluidresulttest.cpp
  src/drt_fluid/xfluidfluidresulttest.cpp
  src/drt_fluid/fluidresulttest.cpp
  src/drt_fluid/time_integration_scheme.cpp
  src/drt_fluid/turbulence_statistic_manager.cpp
  src/drt_fluid/turbulence_statistics_bfs.cpp
  src/drt_fluid/turbulence_statistics_ccy.cpp
  src/drt_fluid/turbulence_statistics_cha.cpp
  src/drt_fluid/turbulence_statistics_ldc.cpp
  src/drt_fluid/turbulence_statistics_mean_general.cpp
  src/drt_fluid/turbulence_statistics_sqc.cpp
  src/drt_fluid/xfluidimplicitintegration.cpp
  src/drt_fluid/xfluidresulttest.cpp
  src/drt_fluid/xfluid.cpp
  src/drt_fluid/xfluidfluid.cpp
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
  src/drt_mlmc/randomfield.cpp 
)
set(OBJS_DRT_MAT
  src/drt_mat/micromaterial_evaluate.cpp
  src/drt_mat/micromaterialgp_static.cpp
)
set(OBJS_BACI_MAIN
  src/drt_lib/drt_init_control.cpp
  src/drt_lib/global_inp_control2.cpp
  src/global_full/global_cal_control.cpp
  src/global_full/global_control.cpp
  src/global_full/baci.cpp
  src/drt_stk/stk_adaptive_main.cpp
)
set(OBJS_DRT_B3II_LIB
  src/drt_beam3ii/beam3ii.cpp
  src/drt_beam3ii/beam3ii_evaluate.cpp
  src/drt_beam3ii/beam3ii_input.cpp
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
  src/drt_structure/str_aux.cpp
  src/drt_structure/strtimada.cpp
  src/drt_structure/strtimada_create.cpp
  src/drt_structure/strtimada_zienxie.cpp
  src/drt_structure/strtimint.cpp
  src/drt_structure/strtimint_ab2.cpp
  src/drt_structure/strtimint_create.cpp
  src/drt_structure/strtimint_expl.cpp
  src/drt_structure/strtimint_gemm.cpp
  src/drt_structure/strtimint_genalpha.cpp
  src/drt_structure/strtimint_impl.cpp
  src/drt_structure/strtimint_impl_nox.cpp
  src/drt_structure/strtimint_noxgroup.cpp
  src/drt_structure/strtimint_noxlinsys.cpp
  src/drt_structure/strtimint_ost.cpp
  src/drt_structure/strtimint_statics.cpp
  src/drt_structure/stru_dyn_nln_drt.cpp
  src/drt_structure/stru_resulttest.cpp
  src/drt_structure/stru_static_drt.cpp
  src/drt_structure/structure_utils_mapextractor.cpp
  src/drt_structure/strudyn_direct.cpp
  src/drt_structure/strugenalpha.cpp
  src/drt_structure/strugenalpha_lsnewton.cpp
)
set(OBJS_POST_DRT_S8CONVERT
  src/post_s8convert/post_drt_s8convert.cpp
)
set(OBJS_LINALG
  src/linalg/aggregation_method.cpp
  src/linalg/aggregation_method_ml.cpp
  src/linalg/aggregation_method_uncoupled.cpp
  src/linalg/amgpreconditioner.cpp
  src/linalg/braesssarazin_smoother.cpp
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
  src/linalg/linalg_solver.cpp
  src/linalg/linalg_sparsematrix.cpp
  src/linalg/linalg_sparsematrixbase.cpp
  src/linalg/linalg_utils.cpp
  src/linalg/saddlepointpreconditioner.cpp
  src/linalg/simpler_operator.cpp
  src/linalg/simpler_operator_ex.cpp
  src/linalg/bgs2x2_operator.cpp
  src/linalg/smoother.cpp
  src/linalg/transfer_operator.cpp
  src/linalg/transfer_operator_pgamg.cpp
  src/linalg/transfer_operator_pgamg2.cpp
  src/linalg/transfer_operator_saamg.cpp
  src/linalg/transfer_operator_tentative.cpp
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
set(OBJS_DRT_B3_LIB
  src/drt_beam3/beam3.cpp
  src/drt_beam3/beam3_evaluate.cpp
  src/drt_beam3/beam3_input.cpp
)
set(OBJS_DRT_F3_LIB
  src/drt_f3/fluid3.cpp
  src/drt_f3/fluid3_boundary_impl.cpp
  src/drt_f3/fluid3_element_boundary.cpp
  src/drt_f3/fluid3_element_boundary_evaluate.cpp
  src/drt_f3/fluid3_evaluate.cpp
  src/drt_f3/fluid3_genalpha_resVMM.cpp
  src/drt_f3/fluid3_input.cpp
  src/drt_f3/fluid3_nurbs.cpp
  src/drt_f3/fluid3_stabilization.cpp
  src/drt_f3/fluid3_weak_dbc.cpp
  src/drt_f3/xfluid3.cpp
  src/drt_f3/xfluid3_evaluate.cpp
  src/drt_f3/xfluid3_input.cpp
  src/drt_f3/xfluid3_line.cpp
  src/drt_f3/xfluid3_line_evaluate.cpp
  src/drt_f3/xfluid3_surface.cpp
  src/drt_f3/xfluid3_surface_evaluate.cpp
  src/drt_f3/xfluid3_sysmat_projection.cpp
  src/drt_f3/xfluid3_sysmat_sigma.cpp
  src/drt_f3/xfluid3_sysmat_tau_pressure.cpp
  src/drt_fluid/time_integration_element.cpp
  src/drt_xdiff3/xdiff3.cpp
  src/drt_xdiff3/xdiff3_evaluate.cpp
  src/drt_xdiff3/xdiff3_input.cpp
  src/drt_xdiff3/xdiff3_line.cpp
  src/drt_xdiff3/xdiff3_line_evaluate.cpp
  src/drt_xdiff3/xdiff3_surface.cpp
  src/drt_xdiff3/xdiff3_surface_evaluate.cpp
  src/drt_xdiff3/xdiff3_sysmat4.cpp
)
set(OBJS_DRT_F3_IMPL_LIB
  src/drt_f3_impl/fluid3_impl.cpp
  src/drt_f3_impl/fluid3_impl_parameter.cpp
)
set(OBJS_DRT_T3_LIB
  src/drt_truss3/truss3.cpp
  src/drt_truss3/truss3_evaluate.cpp
  src/drt_truss3/truss3_input.cpp
)
set(OBJS_DRT_TLM_LIB
  src/drt_trusslm/trusslm.cpp
  src/drt_trusslm/trusslm_evaluate.cpp
  src/drt_trusslm/trusslm_input.cpp
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
  src/drt_contact/strugenalpha_cmt.cpp
)
set(OBJS_STK_LIB
  src/stk_lib/stk_discret.cpp
  src/stk_lib/stk_discret_state.cpp
  src/stk_lib/stk_dofset.cpp
  src/stk_lib/stk_field_comm.cpp
  src/stk_lib/stk_gnuplot.cpp
  src/stk_lib/stk_iterator.cpp
  src/stk_lib/stk_linearsystem.cpp
)
set(OBJS_STK_REFINE
  src/stk_refine/stk_comm.cpp
  src/stk_refine/stk_io.cpp
  src/stk_refine/stk_mesh.cpp
  src/stk_refine/stk_refine.cpp
  src/stk_refine/stk_topology.cpp
  src/stk_refine/stk_unrefine.cpp
  src/stk_refine/stk_utils.cpp
)
set(OBJS_STK_FLUID
  src/stk_fluid/fluid_errorestimate.cpp
  src/stk_fluid/fluid_implicit.cpp
  src/stk_fluid/fluid_problem.cpp
  src/stk_fluid/fluid_resulttest.cpp
)
set(OBJS_STK_STRUCTURE
  src/stk_structure/str_genalpha.cpp
  src/stk_structure/str_problem.cpp
  src/stk_structure/str_resulttest.cpp
)
set(OBJS_STK_FSI
  src/stk_fsi/fsi_problem.cpp
)
set(OBJS_DRT_STK
  src/drt_stk/stk_adaptive.cpp
  src/drt_stk/stk_assemblestrategy.cpp
  src/drt_stk/stk_errorestimate.cpp
  src/drt_stk/stk_fluid.cpp
  src/drt_stk/stk_fluidresulttest.cpp
  src/drt_stk/stk_gmsh.cpp
)
set(OBJS_DRT_CUT
  src/drt_cut/cut_boundarycell.cpp
  src/drt_cut/cut_boundingbox.cpp
  src/drt_cut/cut_creator.cpp
  src/drt_cut/cut_cycle.cpp
  src/drt_cut/cut_edge.cpp
  src/drt_cut/cut_element.cpp
  src/drt_cut/cut_elementhandle.cpp
  src/drt_cut/cut_facet.cpp
  src/drt_cut/cut_facetgraph.cpp
  src/drt_cut/cut_find_cycles.cpp
  src/drt_cut/cut_coloredgraph.cpp
  src/drt_cut/cut_integrationcell.cpp
  src/drt_cut/cut_integrationcellcreator.cpp
  src/drt_cut/cut_levelsetintersection.cpp
  src/drt_cut/cut_levelsetside.cpp
  src/drt_cut/cut_line.cpp
  src/drt_cut/cut_mesh.cpp
  src/drt_cut/cut_meshhandle.cpp
  src/drt_cut/cut_meshintersection.cpp
  src/drt_cut/cut_node.cpp
  src/drt_cut/cut_point.cpp
  src/drt_cut/cut_pointgraph.cpp
  src/drt_cut/cut_pointpool.cpp
  src/drt_cut/cut_side.cpp
  src/drt_cut/cut_sidehandle.cpp
  src/drt_cut/cut_tetmesh.cpp
  src/drt_cut/cut_tetmeshintersection.cpp
  src/drt_cut/cut_volumecell.cpp
)
set(OBJS_CUT_TEST
  tests/cut_test/cut_test_double.cpp
  tests/cut_test/cut_test_facets.cpp
  tests/cut_test/cut_test_graph.cpp
  tests/cut_test/cut_test_generator.cpp
  tests/cut_test/cut_test_intersection.cpp
  tests/cut_test/cut_test_levelset.cpp
  tests/cut_test/cut_test_main.cpp
  tests/cut_test/cut_test_simple.cpp
  tests/cut_test/cut_test_units.cpp
  tests/cut_test/cut_test_utils.cpp
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
  tests/cut_test/cut_test_axel_8.cpp
  tests/cut_test/cut_test_axel_9.cpp
  tests/cut_test/cut_test_axel_10.cpp
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
  src/drt_stru_ale/stru_ale_utils.cpp
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

option (D_WALL1       "wall1" ON)
if (D_WALL1)
  add_definitions("-DD_WALL1")
endif(D_WALL1)

option (D_THERMO      "thermal elements" ON)
if (D_THERMO)
  add_definitions("-DD_THERMO")
endif(D_THERMO)

option (D_SOLID3      "spatial solid elements" ON)
if (D_SOLID3)
  add_definitions("-DD_SOLID3")
endif(D_SOLID3)

option (D_BEAM2       "2D Timoshenko beam element" ON)
if (D_BEAM2)
  add_definitions("-DD_BEAM2")
endif(D_BEAM2)

option (D_BEAM2R      "2D Reissner beam element" ON)
if (D_BEAM2R)
  add_definitions("-DD_BEAM2R")
endif(D_BEAM2R)

option (D_BEAM3       "3D Reissner beam element" ON)
if (D_BEAM3)
  add_definitions("-DD_BEAM3")
endif(D_BEAM3)

option (D_BEAM3II     "3D Reissner beam element of type II" ON)
if (D_BEAM3II)
  add_definitions("-DD_BEAM3II")
endif(D_BEAM3II)

option (D_SMOOTHROD      "3D smooth rod element" ON)
if (D_SMOOTHROD)
  add_definitions("-DD_SMOOTHROD")
endif(D_SMOOTHROD)

option (D_TRUSS3      "3D truss element" ON)
if (D_TRUSS3)
  add_definitions("-DD_TRUSS3")
endif(D_TRUSS3)

option (D_TRUSS2      "2D truss element" ON)
if (D_TRUSS2)
  add_definitions("-DD_TRUSS2")
endif(D_TRUSS2)

option (D_TORSION3    "3D torsion spring element" ON)
if (D_TORSION3)
  add_definitions("-DD_TORSION3")
endif(D_TORSION3)

option (D_TORSION2    "2D torsion spring element" ON)
if (D_TORSION2)
  add_definitions("-DD_TORSION2")
endif(D_TORSION2)

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

###------------------------------------------------------------------ Binaries

# hide libraries

if (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/src)
endif (NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

# Define the libraries. The link order is important!

# Compile element libraries first. There are files that take time to compile.

add_library(ele_f3 ${OBJS_DRT_F3_LIB})
add_library(ele_f3_impl ${OBJS_DRT_F3_IMPL_LIB})
add_library(ele_ale ${OBJS_DRT_ALE_LIB})
add_library(ele_b2 ${OBJS_DRT_B2_LIB})
add_library(ele_b2r ${OBJS_DRT_B2R_LIB})
add_library(ele_b3 ${OBJS_DRT_B3_LIB})
add_library(ele_b3ii ${OBJS_DRT_B3II_LIB})
add_library(ele_bele3 ${OBJS_DRT_BELE3_LIB})
add_library(ele_combust3 ${OBJS_DRT_COMBUST3_LIB})
add_library(ele_sr ${OBJS_DRT_SR_LIB})
add_library(ele_s8 ${OBJS_DRT_S8_LIB})
add_library(ele_scatra ${OBJS_DRT_SCATRA_LIB})
add_library(ele_so ${OBJS_SOH8_LIB})
add_library(ele_t2 ${OBJS_DRT_T2_LIB})
add_library(ele_t3 ${OBJS_DRT_T3_LIB})
add_library(ele_tlm ${OBJS_DRT_TLM_LIB})
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
  ele_bele3
  ele_combust3
  ele_f3
  ele_f3_impl
  ele_sr
  ele_s8
  ele_scatra
  ele_so
  ele_t2
  ele_t3
  ele_tlm
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

add_library(alg_stk_fsi ${OBJS_STK_FSI})

set(MULTIFIELD_LIBRARIES
  alg_elch
  alg_fsi
  alg_loma
  alg_tsi
  alg_stk_fsi
  alg_stru_ale
  )

# the adapter layer
add_library(drt_adapter ${OBJS_DRT_ADAPT})

# algorithm libraries the adapters depend upon
add_library(alg_ale ${OBJS_DRT_ALEALG})
add_library(alg_artnet ${OBJS_DRT_ARTNET})
add_library(alg_combust ${OBJS_DRT_COMBUST})
add_library(alg_contact ${OBJS_DRT_CONTACT} ${OBJS_DRT_BEAMCONTACT})
add_library(alg_fluid ${OBJS_DRT_FLUID})
add_library(alg_invanalysis ${OBJS_INV_ANALYSIS})
add_library(alg_mlmc ${OBJS_MLMC})
add_library(alg_redairways ${OBJS_DRT_RED_AIRWAYS})
add_library(alg_scatra ${OBJS_DRT_SCATRA})
add_library(alg_statmech ${OBJS_DRT_STATMECH})
add_library(alg_structure ${OBJS_DRT_STRUCTURE})
add_library(alg_thermo ${OBJS_DRT_THERMO})

add_library(alg_stk_structure ${OBJS_STK_STRUCTURE})
add_library(alg_stk_fluid ${OBJS_STK_FLUID})

set(SINGLEFIELD_LIBRARIES
  alg_ale
  alg_artnet
  alg_combust
  alg_contact
  alg_fluid
  alg_redairways
  alg_scatra
  alg_statmech
  alg_structure
  alg_invanalysis
  alg_mlmc
  alg_thermo
  alg_stk_structure
  alg_stk_fluid
  )

# general framework
add_library(drt_nurbs ${OBJS_DRT_NURBS_DISCRET})
add_library(drt_constraint ${OBJS_DRT_CONSTRAINT})
add_library(drt_mortar ${OBJS_DRT_MORTAR})

# material definitions
add_library(drt_mat ${OBJS_DRT_MAT_LIB})
add_library(drt_micromat ${OBJS_DRT_MAT} ${OBJS_DRT_STRU_MULTI})

# special library that demands all parobject types
add_library(drt_lib_register ${OBJS_DRT_REGISTER})

# libraries that are used by baci and by filters
add_library(drt_lib ${OBJS_DRT_LIB})
add_library(drt_lib_inpar ${OBJS_DRT_INPAR})
add_library(drt_lib_constraint ${OBJS_DRT_CONSTRAINT_LIB})
add_library(drt_lib_contact ${OBJS_DRT_CONTACT_LIB} ${OBJS_DRT_BEAMCONTACT_LIB})
add_library(drt_lib_intersection ${OBJS_INTERSECTION})
add_library(drt_lib_io ${OBJS_DRT_IO_LIB})
#add_library(drt_lib_fluid ${OBJS_DRT_FLUID_LIB})
add_library(drt_lib_mortar ${OBJS_DRT_MORTAR_LIB})
add_library(drt_lib_nurbs ${OBJS_DRT_NURBS_DISCRET_LIB})
add_library(drt_lib_potential ${OBJS_DRT_POTENTIAL})
add_library(drt_lib_statmech ${OBJS_DRT_STATMECH_LIB})
add_library(drt_lib_surfstress ${OBJS_DRT_SURFSTRESS})
add_library(drt_lib_utils ${OBJS_DRT_UTILS})
add_library(drt_lib_artnet ${OBJS_DRT_ARTNET_LIB})
add_library(drt_lib_redairways ${OBJS_DRT_RED_AIRWAYS_LIB})
add_library(drt_lib_cut ${OBJS_DRT_CUT})

add_library(stk_refine ${OBJS_STK_REFINE})
add_library(stk_drt_lib ${OBJS_STK_LIB})
add_library(drt_lib_stk ${OBJS_DRT_STK})

set(SUPPORT_LIBRARIES
  drt_lib
  drt_lib_artnet 
  drt_lib_redairways
  drt_lib_mortar
  drt_lib_statmech
  drt_lib_constraint
  drt_lib_contact
  # STK based libraries
  drt_lib_stk
  stk_drt_lib
  stk_refine
  # Libraries needed by elements
  drt_lib_surfstress
  drt_lib_potential
  drt_lib_intersection
  drt_lib_cut
  drt_lib_inpar
  drt_lib_io
  drt_lib_nurbs 
  )

# filter framework
add_library(post_drt_common ${OBJS_POST_DRT_COMMON})
add_library(filter_common ${OBJS_FILTER_COMMON})

# linear algebra objects
add_library(linalg ${OBJS_LINALG})

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
  drt_mortar 
  drt_constraint
  drt_lib_register
  ${ELEMENT_LIBRARIES}
  ${SUPPORT_LIBRARIES}
  drt_micromat drt_mat
  drt_lib_utils
  linalg
  ${LIBRARIES})

set(FILTER_LIBRARIES 
  post_drt_common
  drt_lib_register
  ${ELEMENT_LIBRARIES}
  ${SUPPORT_LIBRARIES}
  drt_mat
  drt_lib_utils
  linalg
  filter_common
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

add_executable(cut_test ${OBJS_CUT_TEST})
target_link_libraries(cut_test drt_lib_cut ${QHULL_LIBRARIES} shards epetra teuchos ${LAPACK_LIBRARIES})
add_test(NAME run_cut_test COMMAND $<TARGET_FILE:cut_test>)

add_executable(graph_test ${OBJS_GRAPH_TEST})
target_link_libraries(graph_test ${Boost_LIBRARIES})

if (STKADAPTIVE)

  # Test programs based on stk

  add_library(stk_test_utils tests/stk_test/stk_test_utils.cpp)

  add_executable(stk_test_change_owner tests/stk_test/stk_test_change_owner.cpp)
  target_link_libraries(stk_test_change_owner stk_test_utils ${LIBRARIES})
  add_test(NAME run_stk_test_change_owner COMMAND $<TARGET_FILE:stk_test_change_owner>)

  add_executable(stk_test_corner1 tests/stk_test/stk_test_corner1.cpp)
  target_link_libraries(stk_test_corner1 stk_test_utils ${LIBRARIES})
  add_test(NAME run_stk_test_corner1 COMMAND $<TARGET_FILE:stk_test_corner1>)

  add_executable(stk_test_corner2 tests/stk_test/stk_test_corner2.cpp)
  target_link_libraries(stk_test_corner2 stk_test_utils ${LIBRARIES})
  add_test(NAME run_stk_test_corner2 COMMAND $<TARGET_FILE:stk_test_corner2>)

  add_executable(stk_test_ghosttoshared tests/stk_test/stk_test_ghosttoshared.cpp)
  target_link_libraries(stk_test_ghosttoshared stk_test_utils ${LIBRARIES})
  add_test(NAME run_stk_test_ghosttoshared COMMAND $<TARGET_FILE:stk_test_ghosttoshared>)

  add_executable(stk_test_ownership tests/stk_test/stk_test_ownership.cpp)
  target_link_libraries(stk_test_ownership stk_test_utils ${LIBRARIES})
  add_test(NAME run_stk_test_ownership COMMAND $<TARGET_FILE:stk_test_ownership>)

  add_executable(stk_test_ownership2 tests/stk_test/stk_test_ownership2.cpp)
  target_link_libraries(stk_test_ownership2 stk_test_utils ${LIBRARIES})
  add_test(NAME run_stk_test_ownership2 COMMAND $<TARGET_FILE:stk_test_ownership2>)

endif (STKADAPTIVE)

###------------------------------------------------------------------ Test definitions

macro (baci_test arg nproc restart)
  add_test(NAME ${arg}-p${nproc}
    #CONFIGURATIONS Debug
    COMMAND mpirun -np ${nproc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${arg}.dat xxx)
  #set_tests_properties (${arg}-${nproc} 
  #  PROPERTIES PASS_REGULAR_EXPRESSION "processor 0 finished normally")

  if (${restart})
    add_test(NAME ${arg}-p${nproc}-restart
      #CONFIGURATIONS Debug
      COMMAND mpirun -np ${nproc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${arg}.dat xxx restart=${restart})
    #set_tests_properties (${arg}-${nproc}-r
    #  PROPERTIES PASS_REGULAR_EXPRESSION "processor 0 finished normally")
  endif (${restart})
endmacro (baci_test)

###------------------------------------------------------------------ List of tests
baci_test(artaneurysm_hex 1 "")
baci_test(artaneurysm_hex 2 "")
baci_test(artaneurysm_tet 1 "")
baci_test(artaneurysm_tet 2 "")
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
baci_test(beam3iicontact_augmented 1 "")
baci_test(beam3iicontact_augmented 2 "")
baci_test(beam3iicontact_penalty 1 "")
baci_test(beam3iicontact_penalty 2 "")
baci_test(beam3contactstatmech_augmented 1 "")
baci_test(beam3contactstatmech_augmented 2 "")
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
baci_test(combust_constpressurejump 1 "")
baci_test(combust_constpressurejump 2 "")
baci_test(combust_levelsetcut_horizontal 1 "")
baci_test(combust_levelsetcut_horizontal 2 "")
baci_test(combust_levelsetcut_vertical 1 "")
baci_test(combust_levelsetcut_vertical 2 "")
baci_test(condif1_straight_convection 1 "")
baci_test(condif1_straight_convection 2 "")
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
baci_test(contact2D_hertz_linstatic 1 "")
baci_test(contact2D_hertz_linstatic 2 "")
baci_test(contact2D_lin_augmented 1 "")
baci_test(contact2D_lin_augmented 2 "")
baci_test(contact2D_lin_cmtgena 1 "")
baci_test(contact2D_lin_cmtgena 2 "")
baci_test(contact2D_lin_duallagr 1 "")
baci_test(contact2D_lin_duallagr 2 "")
baci_test(contact2D_lin_penalty 1 "")
baci_test(contact2D_lin_penalty 2 "")
baci_test(contact2D_lin_stdlagr 1 "")
baci_test(contact2D_lin_stdlagr 2 "")
baci_test(contact2D_multibody 1 "")
baci_test(contact2D_multibody 2 "")
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
baci_test(contact3D_lin_cmtgena 1 "")
baci_test(contact3D_lin_cmtgena 2 "")
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
baci_test(contact3D_pwlin_penalty_hex27 1 "")
baci_test(contact3D_pwlin_penalty_hex27 2 "")
baci_test(contact3D_pwlin_penalty_tet10 1 "")
baci_test(contact3D_pwlin_penalty_tet10 2 "")
baci_test(contact3D_quad_penalty_hex27 1 "")
baci_test(contact3D_quad_penalty_hex27 2 "")
baci_test(contact3D_quad_tet10 1 "")
baci_test(contact3D_quad_tet10 2 "")
baci_test(contact3D_slidingblock_augmented 1 "")
baci_test(contact3D_slidingblock_augmented 2 "")
baci_test(contact3D_slidingblock_augmented_tet10 1 "")
baci_test(contact3D_slidingblock_augmented_tet10 2 "")
baci_test(contact3D_slidingblock_lagrange_tet10 1 "")
baci_test(contact3D_slidingblock_lagrange_tet10 2 "")
baci_test(contact3D_slidingblock_lagrange_hex20 1 "")
baci_test(contact3D_slidingblock_lagrange_hex20 2 "")
baci_test(contact3D_slidingblock_cmtgena 1 "")
baci_test(contact3D_slidingblock_cmtgena 2 "")
baci_test(contact3D_slidingblock_coulomb 1 "")
baci_test(contact3D_slidingblock_coulomb 2 "")
baci_test(contact3D_slidingblock_penalty 1 "")
baci_test(contact3D_slidingblock_penalty 2 "")
baci_test(contact3D_wear 1 "")
baci_test(contact3D_wear 2 "")
baci_test(disp_control 1 "")
baci_test(disp_control 2 "")
baci_test(elch_1D_10ele_3ions_stab 1 "")
baci_test(elch_1D_10ele_3ions_stab 2 "")
baci_test(elch_Kwok_Wu 1 "")
baci_test(elch_Kwok_Wu 2 "")
baci_test(elch_Kwok_Wu_BDF2 2 "")
baci_test(elch_RCE_rotsymm_pbc 1 "")
baci_test(elch_RCE_rotsymm_pbc 2 "")
baci_test(elch_natconv 1 2)
baci_test(elch_natconv 2 2)
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
baci_test(elch_pulseplating_hex27_moving_boundary 2 "")
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
baci_test(f2_drivencavity32x32_Re1000_stat 1 "")
baci_test(f2_drivencavity32x32_Re1000_stat 2 "")
baci_test(f2_drivencavity40x40_braess_sarazin_drt 1 "")
baci_test(f2_drivencavity40x40_braess_sarazin_drt 2 "")
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
baci_test(f3_drt_genalpha_initialfield_volforce_pbc 1 "")
baci_test(f3_drt_genalpha_initialfield_volforce_pbc 2 "")
baci_test(f3_drt_impl_initialfield_volforce_pbc 1 "")
baci_test(f3_drt_impl_initialfield_volforce_pbc 2 "")
baci_test(f3_nurbs27_AfGenAlpha 2 "")
baci_test(f3_ew_Q1_27x_drt 1 "")
baci_test(f3_ew_Q1_27x_drt 2 "")
baci_test(f3_onewaycoupled_afgenalpha 2 8)
baci_test(f3_onewaycoupled_genalpha 2 8)
baci_test(f3_onewaycoupled_ost 2 8)
baci_test(f3_simpler_dgs 1 "")
baci_test(f3_simpler_dgs 2 "")
baci_test(fs-kugel-tet 1 "")
baci_test(fs-kugel-tet 2 "")
baci_test(fs-quadrat 1 "")
baci_test(fs-quadrat 2 "")
baci_test(fs-wuerfel-tet 1 "")
baci_test(fs-wuerfel-tet 2 "")
baci_test(fsi-pw-mfsi-be 1 3)
baci_test(fsi-pw-mfsi-be 2 3)
baci_test(fsi-pw-mfsi-be-structsplit 1 3)
baci_test(fsi-pw-mfsi-be-structsplit 2 3)
baci_test(fsi-pw-mfsi-fsiamg-be-structsplit 1 3)
baci_test(fsi-pw-mfsi-fsiamg-be-structsplit 2 3)
baci_test(fsi-pw-mfsi-gena 1 3)
baci_test(fsi-pw-mfsi-gena 2 3)
baci_test(fsi-pw-mfsi-gena-contact 1 "")
baci_test(fsi-pw-mfsi-gena-contact 2 "")
baci_test(fsi_oval3D_sliding_up 1 "")
baci_test(fsi_oval3D_sliding_up 2 "")
baci_test(fsi_ow32x32_drt 1 "")
baci_test(fsi_ow32x32_drt 2 "")
baci_test(fsi_ow32x32_drt_ca 1 "")
baci_test(fsi_ow32x32_drt_ca 2 "")
baci_test(fsi_ow32x32_drt_force 1 "")
baci_test(fsi_ow32x32_drt_force 2 "")
baci_test(fsi_ow32x32_drt_ga_mfsi 1 "")
baci_test(fsi_ow32x32_drt_ga_mfsi 2 "")
baci_test(fsi_ow32x32_drt_gena 1 "")
baci_test(fsi_ow32x32_drt_gena 2 "")
baci_test(fsi_ow32x32_drt_mfsi 1 "")
baci_test(fsi_ow32x32_drt_mfsi 2 "")
baci_test(fsi_ow32x32_sd_drt 1 1)
baci_test(fsi_ow32x32_sd_drt 2 1)
baci_test(fsi_ow32x32_sd_sti_drt 1 1)
baci_test(fsi_ow32x32_sd_sti_drt 2 1)
baci_test(fsi_ow32x32_sti_drt 1 "")
baci_test(fsi_ow32x32_sti_drt 2 "")
baci_test(fsi_ow32x32_sti_drt_mfsi 1 "")
baci_test(fsi_ow32x32_sti_drt_mfsi 2 "")
baci_test(fsi_ow32x32_sti_mtr_drt 1 "")
baci_test(fsi_ow32x32_sti_mtr_drt 2 "")
baci_test(fsi_ow3D_mtr_drt 1 "")
baci_test(fsi_ow3D_mtr_drt 2 "")
baci_test(fsi_potential_LJ 1 "")
baci_test(fsi_potential_LJ_par 2 "")
baci_test(fsi_ring2D_sliding_tot 1 "")
baci_test(fsi_ring2D_sliding_tot 2 "")
baci_test(fsi_ring2D_sliding_up 1 "")
baci_test(fsi_ring2D_sliding_up 2 "")
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
baci_test(mat_balzaniso_plate 1 "")
baci_test(mat_balzaniso_plate 2 "")
baci_test(mat_contchain_remo 1 "")
baci_test(mat_growth_ip 1 " ")
baci_test(mat_growth_ip 2 " ")
baci_test(mat_holzapfelcardio 1 "")
baci_test(mat_holzapfelcardio 2 "")
baci_test(mat_holzharit_remo 1 "")
baci_test(mat_holzharit_remo 2 "")
baci_test(mat_humphreycardio 1 "")
baci_test(mat_humphreycardio 2 "")
baci_test(mat_itskov 1 "")
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
baci_test(mfsi_ow32x32_sti_drt_mtrss 1 "")
baci_test(mfsi_ow32x32_sti_drt_mtrss 2 "")
baci_test(mfsi_ow32x32_sti_drt_ss 1 "")
baci_test(mfsi_ow32x32_sti_drt_ss 2 "")
baci_test(mfsi_ow32x32_sti_drt_mtrfs 1 "")
baci_test(mfsi_ow32x32_sti_drt_mtrfs 2 "")
baci_test(mfsi_ow32x32_sti_drt_fs 1 "")
baci_test(mfsi_ow32x32_sti_drt_fs 2 "")
baci_test(mfsi_ow3D_mtr_drt 1 "")
baci_test(mfsi_ow3D_mtr_drt 2 "")
baci_test(mfsi_oval3D_sliding_up_fs 1 "")
baci_test(mfsi_oval3D_sliding_up_fs 2 "")
baci_test(mfsi_ves3D_slid_curr_fs 1 "")
baci_test(mfsi_ves3D_slid_curr_fs 2 "")
baci_test(mfsi_ves3D_slid_curr_ss 1 "")
baci_test(mfsi_ves3D_slid_curr_ss 2 "")
baci_test(mfsi_ves3D_slid_rot_fs 1 "")
baci_test(mfsi_ves3D_slid_rot_fs 2 "")
baci_test(mfsi_ves3D_slid_rot_ss 1 "")
baci_test(mfsi_ves3D_slid_rot_ss 2 "")
baci_test(orthopressure_shell8_eas_ans_sphere 1 8)
baci_test(orthopressure_shell8_eas_ans_sphere 2 8)
baci_test(orthopressure_solidh8_eas_sphere 1 "")
baci_test(orthopressure_solidh8_eas_sphere 2 "")
baci_test(orthopressure_solidh8_sphere 1 "")
baci_test(orthopressure_solidh8_sphere 2 "")
baci_test(orthopressure_solidw6_sphere 1 "")
baci_test(orthopressure_solidw6_sphere 2 "")
baci_test(potential_surf_pointapprox_VDW_3D 1 "")
baci_test(potential_surf_pointapprox_VDW_3D 2 "")
baci_test(potential_vol_LJ_2D 1 "")
baci_test(potential_vol_LJ_2D 2 "")
baci_test(potential_vol_LJ_3D 1 "")
baci_test(potential_vol_LJ_3D 2 "")
baci_test(prestress_fsi 1 "")
baci_test(prestress_fsi 2 "")
baci_test(prestress_fsi_id 1 "")
baci_test(prestress_fsi_id 2 "")
baci_test(prestress_hex8 1 "")
baci_test(prestress_hex8 2 "")
baci_test(prestress_tet4 1 "")
baci_test(prestress_tet4 2 "")
baci_test(prestress_weg6 1 "")
baci_test(prestress_weg6 2 "")
baci_test(scatra_gaussian_hill_pbc 1 "")
baci_test(scatra_gaussian_hill_pbc 2 "")
baci_test(scatra_gaussian_hill_pbc_gena_cons 1 "")
baci_test(scatra_gaussian_hill_pbc_gena_cons 2 "")
baci_test(scatra_gaussian_hill_pbc_genalpha 1 "")
baci_test(scatra_gaussian_hill_pbc_genalpha 2 "")
baci_test(scatra_graetz_problem_system 1 5)
baci_test(scatra_graetz_problem_system 2 5)
baci_test(scatra_graetz_problem_system_ost 2 5)
baci_test(scatra_nurbs9_channel 1 "")
baci_test(scatra_nurbs9_channel 2 5)
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
# baci_test(sohex8fbar_necking_plastic 1 "")
# baci_test(sohex8fbar_necking_plastic 2 "")
baci_test(sohex8p1j1_cooks_nln_drt 1 "")
baci_test(sohex8p1j1_cooks_nln_drt 2 "")
# baci_test(sohex8p1j1_necking_plastic 1 "")
# baci_test(sohex8p1j1_necking_plastic 2 "")
baci_test(solidh8_mat_neohooke_drt 1 "")
baci_test(solidh8_mat_neohooke_drt 2 "")
baci_test(sosh8_freeflying_ruler 1 "")
baci_test(sosh8_freeflying_ruler 2 "")
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
baci_test(sotet10_cooks_nl 1 "")
baci_test(sotet10_cooks_nl 2 "")
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
baci_test(strufem_ptet4_disp_cantilever 1 "")
baci_test(strufem_ptet4_disp_cantilever 2 "")
baci_test(strufem_shell8_eas_cantilever 1 "")
baci_test(strufem_shell8_eas_cantilever 2 "")
baci_test(strufem_sosh8_eas_cantilever 1 "")
baci_test(strufem_sosh8_eas_cantilever 2 "")
baci_test(strufem_tet10_disp_cantilever 1 "")
baci_test(strufem_tet10_disp_cantilever 2 "")
baci_test(strufem_tet4_disp_cantilever 1 "")
baci_test(strufem_tet4_disp_cantilever 2 "")
baci_test(strufem_wedge6_disp_cantilever 1 "")
baci_test(strufem_wedge6_disp_cantilever 2 "")
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
baci_test(tsi_contact3D_twobodies 1 "")
baci_test(tsi_contact3D_twobodies 2 "")
baci_test(tsi_contact3D_heatproduction 1 "")
baci_test(tsi_contact3D_heatproduction 2 "")
baci_test(tsi_contact3D_twowaycoupling 1 "")
baci_test(tsi_contact3D_twowaycoupling 2 "")
baci_test(tsi_heatconvection 1 "")
baci_test(tsi_heatconvection 2 "")
baci_test(tsi_lincompression_iterstaggdisp 1 "")
baci_test(tsi_lincompression_iterstaggtemp 1 "")
baci_test(tsi_lincompression_sequstaggdisp 1 "")
baci_test(tsi_lincompression_sequstaggtemp 1 "")
baci_test(tsi_lincompression_1waydisp 1 "")
baci_test(tsi_lindilatation_geolin 1 "")
baci_test(viscoAnisotropic_anisocreep 1 60)
baci_test(viscoAnisotropic_anisocreep 2 60)
baci_test(viscoAnisotropic_creep 1 60)
baci_test(viscoAnisotropic_creep 2 60)
baci_test(viscoNeoHooke_creep 1 60)
baci_test(viscoNeoHooke_creep 2 60)
baci_test(w1_dyn_dirich_drt 1 "")
baci_test(w1_dyn_dirich_drt 2 "")
baci_test(w1q4_cantilever_eas 1 "")
baci_test(w1q4_cantilever_eas 2 "")
baci_test(w1q4_patch_eas 1 "")
baci_test(w1q4_patch_eas 2 "")
baci_test(w1q9_dyn_rot_drt_gemm 1 "")
baci_test(w1q9_dyn_rot_drt_gena 1 "")
baci_test(womersley_f3 2 "")

if (D_BEAM3II)
  baci_test(beam3iitest1linear 1 "")
  baci_test(beam3iitest1linear 2 "")
  baci_test(beam3iitest2linear 1 "")
  baci_test(beam3iitest2linear 2 "")
endif (D_BEAM3II)

baci_test(inversedesign_wedge6 1 "")
baci_test(inversedesign_wedge6 2 "")

baci_test(fluidfluidCoupling_hex20_constinflow 1 "")
#  baci_test(fluidfluidCoupling_hex20_parainflow_ydir 1 "")
#  baci_test(fluidfluidCoupling_hex20_parainflow_ydir_cond 1 "")
baci_test(fluidfluidCoupling_simple 1 "")
baci_test(fxf_shearflow_stat 1 "")
baci_test(fxf_shearflow_stat_cond 1 "")
baci_test(xff_channel_simple_stat 1 "")
baci_test(two_phase_couette_flow 1 "")
baci_test(two_phase_rayleigh_taylor_instability_32x128 1 4)
baci_test(two_phase_stat_bubble_surftens_20x20x20 2 "")
baci_test(two_phase_stat_bubble_surftens_50x50x1 1 "")
baci_test(xf3_beltrami_tet4_drt 1 2)
baci_test(xf3_beltrami_tet4_drt 2 2)
baci_test(xf3_drivencavity_stat_drt 1 "")
baci_test(xf3_drivencavity_stat_drt 2 "")
baci_test(xf3_ew_Q1_27x_drt 1 15)
baci_test(xf3_ew_Q1_27x_drt 2 15)
baci_test(xf3_flow_over_corner_fixpoint 1 "")
baci_test(xf3_flow_over_corner_fixpoint 2 "")
baci_test(xf3_flow_over_corner_newton 1 "")
baci_test(xf3_multibody_hex8 1 "")
baci_test(xf3_multibody_hex8 2 "")
baci_test(xf3_multibody_hex8_tp_condense_no 1 "")
baci_test(xf3_multibody_hex8_tp_condense_no 2 "")
baci_test(xf3_multibody_hex8_tp_condense_yes 1 "")
baci_test(xf3_multibody_hex8_tp_condense_yes 2 "")
baci_test(xf3_multibody_hex8_2 1 "")
baci_test(xf3_pure_shear_forces 1 "")
baci_test(xf3_pure_shear_forces 2 "")
baci_test(xf3_pure_shear_forces_bdf2 1 10)
baci_test(xf3_pure_shear_forces_bdf2 2 10)
baci_test(xf3_pure_shear_forces_ost 1 15)
baci_test(xf3_pure_shear_forces_ost 2 15)
baci_test(xf3_pure_shear_forces_tp 1 "")
baci_test(xf3_pure_shear_forces_tp 2 "")
baci_test(xf3_pure_shear_forces_tp_ost 1 15)
baci_test(xf3_pure_shear_forces_tp_ost 2 15)
baci_test(xf3_pure_shear_forces_tp_ost_no_condense 1 15)
baci_test(xf3_pure_shear_forces_tp_ost_no_condense 2 15)
baci_test(xfsi_RingInChannel3D 1 "")
baci_test(xfsi_RingInChannel3D 2 "")
baci_test(xfsi_bending_2d_stationary_hex20 1 "")
baci_test(xfsi_bending_2d_stationary_hex20 2 "")
baci_test(xfsi_bending_2d_stationary_hex8 1 "")
baci_test(xfsi_bending_2d_stationary_hex8 2 "")
baci_test(xfsi_fluid_pusher_BE 1 15)
baci_test(xfsi_fluid_pusher_bdf2 1 "")
baci_test(xfsi_fluid_pusher_bdf2_monolithic 1 "")
baci_test(xfsi_fluid_pusher_bdf2_monolithic_schmu 1 "")
baci_test(xfsi_fluid_pusher_ost 1 15)
baci_test(xfsi_fluid_pusher_ost_densmod 1 15)
baci_test(xfsi_fluid_pusher_ost_monolithic 1 15)
baci_test(xfsi_fluid_pusher_ost_projection 1 15)
baci_test(xfsi_fluid_pusher_stationary_hex8_monolithic 1 "")
baci_test(xfsi_fluid_pusher_stationary_hex8_partit 1 "")
baci_test(xfsi_fluid_pusher_stationary_hex8_partit 2 "")
baci_test(xfsi_multibody_hex20 1 "")
baci_test(xfsi_multibody_hex20 2 "")
baci_test(xfsi_multibody_hex20_tp_DLM_yes 1 "")
baci_test(xfsi_multibody_hex20_tp_DLM_yes 2 "")
baci_test(xfsi_multibody_hex8_monolithic 1 "")
baci_test(xfsi_shooter_thin 1 "")

# remove any output files from our tests
add_test(NAME test_cleanup COMMAND sh -c "rm -v xxx* core.*")
