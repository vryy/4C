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
set(OBJS_DRT_MAT_LIB
    src/drt_matelast/elast_aniso_structuraltensor_strategy.cpp
    src/drt_mat/aaa_mixedeffects.cpp
    src/drt_mat/aaagasser.cpp
    src/drt_mat/aaaneohooke.cpp
    src/drt_mat/aaaneohooke_stopro.cpp
    src/drt_mat/aaaraghavanvorp_damage.cpp
    src/drt_mat/anisotropy.cpp
    src/drt_mat/anisotropy_cylinder_coordinate_system_manager.cpp
    src/drt_mat/anisotropy_extension.cpp
    src/drt_mat/anisotropy_extension_cylinder_cosy.cpp
    src/drt_mat/anisotropy_extension_default.cpp
    src/drt_mat/anisotropy_utils.cpp
    src/drt_mat/arrhenius_pv.cpp
    src/drt_mat/arrhenius_spec.cpp
    src/drt_mat/arrhenius_temp.cpp
    src/drt_mat/beam_elasthyper.cpp
    src/drt_mat/beam_elasthyper_parameter.cpp
    src/drt_mat/beam3r_plasticity.cpp
    src/drt_mat/carreauyasuda.cpp
    src/drt_mat/cnst_1d_art.cpp
    src/drt_mat/compogden.cpp
    src/drt_mat/constraintmixture.cpp
    src/drt_mat/crosslinkermat.cpp
    src/drt_mat/crystal_plasticity.cpp
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
    src/drt_mat/fluid_weakly_compressible.cpp
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
    src/drt_mat/membrane_elasthyper.cpp
    src/drt_mat/membrane_elasthyper_service.cpp
    src/drt_mat/membrane_active_strain.cpp
    src/drt_mat/growthremodel_elasthyper.cpp
    src/drt_mat/growth.cpp
    src/drt_mat/growth_law.cpp
    src/drt_mat/herschelbulkley.cpp
    src/drt_mat/ion.cpp
    src/drt_mat/inelastic_defgrad_factors.cpp
    src/drt_mat/lin_elast_1D.cpp
    src/drt_mat/lubrication_law.cpp
    src/drt_mat/lubrication_mat.cpp
    src/drt_mat/material.cpp
    src/drt_mat/material_service.cpp
    src/drt_mat/matlist.cpp
    src/drt_mat/matlist_reactions.cpp
    src/drt_mat/matpar_bundle.cpp
    src/drt_mat/matpar_material.cpp
    src/drt_mat/matpar_parameter.cpp
    src/drt_mat/micromaterial.cpp
    src/drt_mat/mixfrac.cpp
    src/drt_mat/modpowerlaw.cpp
    src/drt_mat/multiplicative_split_defgrad_elasthyper.cpp
    src/drt_mat/multiplicative_split_defgrad_elasthyper_service.H
    src/drt_mat/myocard.cpp
    src/drt_mat/myocard_fitzhugh_nagumo.cpp
    src/drt_mat/myocard_inada.cpp
    src/drt_mat/myocard_minimal.cpp
    src/drt_mat/myocard_san_garny.cpp
    src/drt_mat/myocard_tentusscher.cpp
    src/drt_mat/myocard_tools.cpp
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
    src/drt_mat/scatra_multiscale.cpp
    src/drt_mat/scatra_reaction_mat.cpp
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
    src/drt_mat/viscoplastic_no_yield_surface.cpp
    src/drt_mat/yoghurt.cpp
    src/drt_matelast/elast_anisoactivestress_evolution.cpp
    src/drt_matelast/elast_coupanisoexpo.cpp
    src/drt_matelast/elast_coupanisoexpoactive.cpp
    src/drt_matelast/elast_coupanisoexpobase.cpp
    src/drt_matelast/elast_coupanisoexposhear.cpp
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
    src/drt_matelast/elast_isomuscle_blemker.cpp
    src/drt_matelast/elast_isoneohooke.cpp
    src/drt_matelast/elast_isoogden.cpp
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
    src/drt_mat/mixture_elasthyper.cpp
    src/drt_mat/hemoglobin_0d_O2_saturation.cpp
    src/drt_mat/air_0d_O2_saturation.cpp
    src/drt_mat/particle_material_base.cpp
    src/drt_mat/particle_material_thermo.cpp
    src/drt_mat/particle_material_sph_fluid.cpp
    src/drt_mat/particle_material_sph_boundary.cpp
    src/drt_mat/particle_material_dem.cpp
    src/drt_mat/particle_wall_material_dem.cpp
    src/drt_mat/activefiber.cpp
    src/drt_matelast/visco_generalizedgenmax.cpp
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
    src/drt_inpar/drt_validcontactconstitutivelaw.cpp
    src/drt_inpar/drt_validmaterials.cpp
    src/drt_inpar/drt_validparameters.cpp
    src/drt_inpar/inpar_cardiovascular0d.cpp
    src/drt_inpar/inpar_ale.cpp
    src/drt_inpar/inpar_beamcontact.cpp
    src/drt_inpar/inpar_beaminteraction.cpp
    src/drt_inpar/inpar_beampotential.cpp
    src/drt_inpar/inpar_beam_to_solid.cpp
    src/drt_inpar/inpar_browniandyn.cpp
    src/drt_inpar/inpar_geometry_pair.cpp
    src/drt_inpar/inpar_geometric_search.cpp
    src/drt_inpar/inpar_bio.cpp
    src/drt_inpar/inpar_cardiac_monodomain.cpp
    src/drt_inpar/inpar_particle.cpp
    src/drt_inpar/inpar_pasi.cpp
    src/drt_inpar/inpar_contact.cpp
    src/drt_inpar/inpar_cut.cpp
    src/drt_inpar/inpar_ehl.cpp
    src/drt_inpar/inpar_elemag.cpp
    src/drt_inpar/inpar_elch.cpp
    src/drt_inpar/inpar_fbi.cpp
    src/drt_inpar/inpar_fluid.cpp
    src/drt_inpar/inpar_fpsi.cpp
    src/drt_inpar/inpar_fs3i.cpp
    src/drt_inpar/inpar_fsi.cpp
    src/drt_inpar/inpar_immersed.cpp
    src/drt_inpar/inpar_invanalysis.cpp
    src/drt_inpar/inpar_io.cpp
    src/drt_inpar/inpar_IO_monitor_structure_dbc.cpp
    src/drt_inpar/inpar_IO_runtime_vtk_output.cpp
    src/drt_inpar/inpar_IO_runtime_vtk_output_fluid.cpp
    src/drt_inpar/inpar_IO_runtime_vtk_output_structure.cpp
    src/drt_inpar/inpar_IO_runtime_vtk_output_structure_beams.cpp
    src/drt_inpar/inpar_IO_runtime_vtp_output_structure.cpp
    src/drt_inpar/inpar_levelset.cpp
    src/drt_inpar/inpar_loca_continuation.cpp
    src/drt_inpar/inpar_lubrication.cpp
    src/drt_inpar/inpar_binningstrategy.cpp
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
    src/drt_inpar/inpar_ssti.cpp
    src/drt_inpar/inpar_sti.cpp
    src/drt_inpar/inpar_structure.cpp
    src/drt_inpar/inpar_thermo.cpp
    src/drt_inpar/inpar_topopt.cpp
    src/drt_inpar/inpar_tsi.cpp
    src/drt_inpar/inpar_twophase.cpp
    src/drt_inpar/inpar_volmortar.cpp
    src/drt_inpar/inpar_wear.cpp
    src/drt_inpar/inpar_xfem.cpp
    src/drt_tutorial/inpar_tutorial.cpp
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
    src/drt_xfem/xfem_xfluid_contact_communicator.cpp
    )

set(OBJS_POST_DRT_MONITOR src/post_monitor/post_drt_monitor.cpp)

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
    src/drt_so3/so_element_service.cpp
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
set(OBJS_DRT_UTILS
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
    src/drt_lib/drt_utils_parallel.cpp
    src/drt_lib/epetra_utils.cpp
    src/drt_beaminteraction/drt_utils_parallel_proctoproc.cpp
    src/drt_lib/drt_utils_rebalancing.cpp
    src/drt_structure_new/str_functions.cpp
    src/pss_full/pss_am.c
    src/pss_full/pss_table.c
    src/pss_full/pss_table_iter.c
    src/drt_lib/voigt_notation.cpp
    )
set(OBJS_DRT_SCATRA
    src/drt_scatra/scatra_algorithm.cpp
    src/drt_scatra/scatra_dyn.cpp
    src/drt_scatra/scatra_cardiac_monodomain_dyn.cpp
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
    src/drt_scatra/scatra_timint_meshtying_strategy_artery.cpp
    src/drt_scatra/scatra_timint_elch.cpp
    src/drt_scatra/scatra_timint_elch_scheme.cpp
    src/drt_scatra/scatra_timint_elch_service.cpp
    src/drt_scatra/scatra_timint_poromulti.cpp
    src/drt_scatra/scatra_utils_clonestrategy.cpp
    src/drt_scatra/scatra_utils.cpp
    src/drt_scatra/turbulence_hit_initial_scalar_field.cpp
    src/drt_scatra/turbulence_hit_scalar_forcing.cpp
    src/drt_scatra/scatra_timint_hdg.cpp
    src/drt_scatra/scatra_timint_stat_hdg.cpp
    src/drt_scatra/scatra_timint_cardiac_monodomain_scheme_hdg.cpp
    src/drt_scatra/nox_nln_scatra_linearsystem.cpp
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
set(OBJS_DRT_FBI
    src/drt_fbi/fbi_calc_utils.cpp
    src/drt_fbi/constraintenforcer_fbi.cpp
    src/drt_fbi/ad_fbi_constraintbridge_penalty.cpp
    src/drt_fbi/ad_fbi_constraintbridge.cpp
    src/drt_fbi/beam_to_fluid_assembly_manager_factory.cpp
    src/drt_fbi/beam_to_fluid_meshtying_pair_base.cpp
    src/drt_fbi/beam_to_fluid_meshtying_pair_factory.cpp
    src/drt_fbi/beam_to_fluid_meshtying_pair_gauss_point.cpp
    src/drt_fbi/beam_to_fluid_meshtying_pair_mortar.cpp
    src/drt_fbi/beam_to_fluid_meshtying_params.cpp
    src/drt_fbi/beam_to_fluid_meshtying_vtk_output_params.cpp
    src/drt_fbi/beam_to_fluid_meshtying_vtk_output_writer.cpp
    src/drt_fbi/beam_to_fluid_mortar_manager.cpp
    src/drt_fbi/constraintenforcer_fbi_factory.cpp
    src/drt_fbi/constraintenforcer_fbi_penalty.cpp
    src/drt_fbi/fluid_assembly_strategy.cpp
    src/drt_fbi/fluidblockmatrix_assembly_strategy.cpp
    src/drt_fbi/immersed_geometry_coupler_fbi.cpp
    src/drt_fbi/immersed_geometry_coupler_fbi_binning.cpp
    src/drt_fbi/immersedcoupler_factory.cpp
    src/drt_fbi/partitioned_penaltycoupling_assembly_manager.cpp
    src/drt_fbi/partitioned_penaltycoupling_assembly_manager_direct.cpp
    src/drt_fbi/partitioned_penaltycoupling_assembly_manager_indirect.cpp
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
    src/drt_fsi/fsi_dirichletneumannslideale.cpp
    src/drt_fsi/fsi_dirichletneumann_disp.cpp
    src/drt_fsi/fsi_dirichletneumann_vel.cpp
    src/drt_fsi/fsi_dirichletneumann_volcoupl.cpp
    src/drt_fsi/fsi_dyn.cpp
    src/drt_fsi/fsi_fluid_ale.cpp
    src/drt_fsi/fsi_lungmonolithic.cpp
    src/drt_fsi/fsi_lungmonolithic_fluidsplit.cpp
    src/drt_fsi/fsi_lungmonolithic_structuresplit.cpp
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

set(OBJS_DRT_REGISTER src/drt_lib/drt_parobjectregister.cpp)

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
    src/drt_scatra_ele/scatra_ele_boundary_calc_elch_electrode_utils.cpp
    src/drt_scatra_ele/scatra_ele_boundary_calc_elch_electrode_growth.cpp
    src/drt_scatra_ele/scatra_ele_boundary_calc_elch_electrode_growth_utils.cpp
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
    src/drt_scatra_ele/scatra_ele_calc_advanced_reaction.cpp
    src/drt_scatra_ele/scatra_ele_calc_refconc_reac.cpp
    src/drt_scatra_ele/scatra_ele_calc_aniso.cpp
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
    src/drt_scatra_ele/scatra_ele_parameter_boundary.cpp
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
    src/drt_sti/sti_monolithic_evaluate_OffDiag.cpp
    src/drt_sti/sti_partitioned.cpp
    src/drt_sti/sti_resulttest.cpp
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
set(OBJS_DRT_NURBS_DISCRET_LIB
    src/drt_fem_general/drt_utils_bspline.cpp
    src/drt_fem_general/drt_utils_nurbs_shapefunctions.cpp
    src/drt_nurbs_discret/drt_control_point.cpp
    src/drt_nurbs_discret/drt_knotvector.cpp
    src/drt_nurbs_discret/drt_nurbs_discret.cpp
    )

set(OBJS_DRT_TWO_PHASE
    src/drt_two_phase_flow/two_phase_algorithm.cpp src/drt_two_phase_flow/two_phase_dyn.cpp
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
set(OBJS_INV_ANALYSIS src/drt_inv_analysis/invana_cal_drt.cpp)
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
set(OBJS_STATINV_ANALYSIS_LIB src/drt_inv_analysis/particle_data.cpp)
set(OBJS_DRT_MAT src/drt_mat/micromaterial_evaluate.cpp src/drt_mat/micromaterialgp_static.cpp)
set(OBJS_DRT_MAT_SCATRA_MULTI
    src/drt_mat/scatra_multiscale_evaluate.cpp src/drt_mat/scatra_multiscale_gp.cpp
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
set(OBJS_PRE_LOCSYS src/pre_locsys/pre_locsys.cpp)
set(OBJS_DRT_SURFSTRESS src/drt_surfstress/drt_surfstress_manager.cpp)
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
    src/drt_mortar/mortar_interface_utils.cpp
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
    src/drt_structure/stru_dyn_nln_drt.cpp
    src/drt_structure_new/gauss_point_data_output_manager.cpp
    src/drt_structure_new/str_utils.cpp
    src/drt_structure_new/str_dbc.cpp
    src/drt_structure_new/str_monitor_dbc.cpp
    src/drt_structure_new/str_timint_basedataio_monitor_dbc.cpp
    src/drt_structure_new/str_resulttest.cpp
    src/drt_structure_new/str_factory.cpp
    src/drt_structure_new/str_impl_generic.cpp
    src/drt_structure_new/str_impl_statics.cpp
    src/drt_structure_new/str_impl_genalpha.cpp
    src/drt_structure_new/str_impl_genalpha_liegroup.cpp
    src/drt_structure_new/str_impl_gemm.cpp
    src/drt_structure_new/str_impl_ost.cpp
    src/drt_structure_new/str_impl_prestress.cpp
    src/drt_structure_new/str_integrator.cpp
    src/drt_structure_new/str_model_evaluator.cpp
    src/drt_structure_new/str_model_evaluator_data.cpp
    src/drt_structure_new/str_model_evaluator_contactdata.cpp
    src/drt_structure_new/str_model_evaluator_factory.cpp
    src/drt_structure_new/str_model_evaluator_generic.cpp
    src/drt_structure_new/str_model_evaluator_structure.cpp
    src/drt_structure_new/str_model_evaluator_cardiovascular0d.cpp
    src/drt_structure_new/str_model_evaluator_lagpenconstraint.cpp
    src/drt_structure_new/str_model_evaluator_springdashpot.cpp
    src/drt_structure_new/str_model_evaluator_contact.cpp
    src/drt_structure_new/str_model_evaluator_meshtying.cpp
    src/drt_struct_ale/struct_ale_str_model_evaluator.cpp
    src/drt_structure_new/str_model_evaluator_multiphysics.cpp
    src/drt_structure_new/str_nln_solver_factory.cpp
    src/drt_structure_new/str_nln_solver_utils.cpp
    src/drt_structure_new/str_nln_solver_generic.cpp
    src/drt_structure_new/str_nln_solver_fullnewton.cpp
    src/drt_structure_new/str_nln_solver_nox.cpp
    src/drt_structure_new/str_nln_solver_ptc.cpp
    src/drt_structure_new/nox_nln_str_linearsystem.cpp
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
    )
set(OBJS_POST_DRT_S8CONVERT src/post_s8convert/post_drt_s8convert.cpp)
set(OBJS_DRT_PARTICLE_ALGORITHM
    src/drt_particle_algorithm/particle_algorithm_utils.cpp
    src/drt_particle_algorithm/particle_algorithm.cpp
    src/drt_particle_algorithm/particle_dirichlet_bc.cpp
    src/drt_particle_algorithm/particle_gravity.cpp
    src/drt_particle_algorithm/particle_initial_field.cpp
    src/drt_particle_algorithm/particle_input_generator.cpp
    src/drt_particle_algorithm/particle_result_test.cpp
    src/drt_particle_algorithm/particle_sim.cpp
    src/drt_particle_algorithm/particle_temperature_bc.cpp
    src/drt_particle_algorithm/particle_timint.cpp
    src/drt_particle_algorithm/particle_viscous_damping.cpp
    )
set(OBJS_DRT_PARTICLE_ENGINE
    src/drt_particle_engine/particle_communication_utils.cpp
    src/drt_particle_engine/particle_container.cpp
    src/drt_particle_engine/particle_container_bundle.cpp
    src/drt_particle_engine/particle_engine.cpp
    src/drt_particle_engine/particle_runtime_vtp_writer.cpp
    src/drt_particle_engine/particle_unique_global_id.cpp
    )
set(OBJS_DRT_PARTICLE_ENGINE_LIB
    src/drt_particle_engine/particle_enums.cpp src/drt_particle_engine/particle_object.cpp
    )
set(OBJS_DRT_PARTICLE_INTERACTION
    src/drt_particle_interaction/particle_interaction_base.cpp
    src/drt_particle_interaction/particle_interaction_dem_adhesion_law.cpp
    src/drt_particle_interaction/particle_interaction_dem_adhesion_surface_energy.cpp
    src/drt_particle_interaction/particle_interaction_dem_adhesion.cpp
    src/drt_particle_interaction/particle_interaction_dem_contact_normal.cpp
    src/drt_particle_interaction/particle_interaction_dem_contact_rolling.cpp
    src/drt_particle_interaction/particle_interaction_dem_contact_tangential.cpp
    src/drt_particle_interaction/particle_interaction_dem_contact.cpp
    src/drt_particle_interaction/particle_interaction_dem_history_pairs.cpp
    src/drt_particle_interaction/particle_interaction_dem_neighbor_pairs.cpp
    src/drt_particle_interaction/particle_interaction_dem.cpp
    src/drt_particle_interaction/particle_interaction_material_handler.cpp
    src/drt_particle_interaction/particle_interaction_runtime_writer.cpp
    src/drt_particle_interaction/particle_interaction_sph_artificialviscosity.cpp
    src/drt_particle_interaction/particle_interaction_sph_boundary_particle.cpp
    src/drt_particle_interaction/particle_interaction_sph_density_correction.cpp
    src/drt_particle_interaction/particle_interaction_sph_density.cpp
    src/drt_particle_interaction/particle_interaction_sph_equationofstate_bundle.cpp
    src/drt_particle_interaction/particle_interaction_sph_equationofstate.cpp
    src/drt_particle_interaction/particle_interaction_sph_heatloss_evaporation.cpp
    src/drt_particle_interaction/particle_interaction_sph_heatsource.cpp
    src/drt_particle_interaction/particle_interaction_sph_kernel.cpp
    src/drt_particle_interaction/particle_interaction_sph_momentum_formulation.cpp
    src/drt_particle_interaction/particle_interaction_sph_momentum.cpp
    src/drt_particle_interaction/particle_interaction_sph_neighbor_pairs.cpp
    src/drt_particle_interaction/particle_interaction_sph_open_boundary.cpp
    src/drt_particle_interaction/particle_interaction_sph_phase_change.cpp
    src/drt_particle_interaction/particle_interaction_sph_pressure.cpp
    src/drt_particle_interaction/particle_interaction_sph_rigid_particle_contact.cpp
    src/drt_particle_interaction/particle_interaction_sph_surface_tension.cpp
    src/drt_particle_interaction/particle_interaction_sph_surface_tension_barrier_force.cpp
    src/drt_particle_interaction/particle_interaction_sph_surface_tension_interface_viscosity.cpp
    src/drt_particle_interaction/particle_interaction_sph_surface_tension_recoilpressure_evaporation.cpp
    src/drt_particle_interaction/particle_interaction_sph_temperature.cpp
    src/drt_particle_interaction/particle_interaction_sph_virtual_wall_particle.cpp
    src/drt_particle_interaction/particle_interaction_sph.cpp
    )
set(OBJS_DRT_PARTICLE_RIGIDBODY
    src/drt_particle_rigidbody/particle_rigidbody_affiliation_pairs.cpp
    src/drt_particle_rigidbody/particle_rigidbody_datastate.cpp
    src/drt_particle_rigidbody/particle_rigidbody_result_test.cpp
    src/drt_particle_rigidbody/particle_rigidbody_runtime_vtp_writer.cpp
    src/drt_particle_rigidbody/particle_rigidbody.cpp
    )
set(OBJS_DRT_PARTICLE_WALL
    src/drt_particle_wall/particle_wall_datastate.cpp
    src/drt_particle_wall/particle_wall_discretization_runtime_vtu_writer.cpp
    src/drt_particle_wall/particle_wall_result_test.cpp
    src/drt_particle_wall/particle_wall.cpp
    )

set(OBJS_LOCA_CONTINUATION
    src/loca_continuation/loca_nln_problem.cpp
    src/loca_continuation/loca_nln_group.cpp
    src/loca_continuation/loca_nln_stepper.cpp
    src/loca_continuation/loca_nln_statustest_factory.cpp
    )

set(OBJS_DRT_NURBS_DISCRET src/drt_nurbs_discret/drt_apply_nurbs_initial_condition.cpp)

set(OBJS_DRT_T3_LIB
    src/drt_truss3/truss3.cpp
    src/drt_truss3/truss3_evaluate.cpp
    src/drt_truss3/truss3_input.cpp
    src/drt_truss3/truss3_scatra.cpp
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
    src/drt_wear/wear_algorithm.cpp src/drt_wear/wear_dyn.cpp src/drt_wear/wear_partitioned.cpp
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
    src/drt_poromultiphase_scatra/poromultiphase_scatra_artery_coupling_nonconforming.cpp
    src/drt_poromultiphase_scatra/poromultiphase_scatra_artery_coupling_nodebased.cpp
    src/drt_poromultiphase_scatra/poromultiphase_scatra_artery_coupling_nodetopoint.cpp
    src/drt_poromultiphase_scatra/poromultiphase_scatra_artery_coupling_linebased.cpp
    src/drt_poromultiphase_scatra/poromultiphase_scatra_artery_coupling_surfbased.cpp
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
    src/drt_ssi/ssi_manifold_utils.cpp
    src/drt_ssi/ssi_monolithic.cpp
    src/drt_ssi/ssi_monolithic_assemble_strategy.cpp
    src/drt_ssi/ssi_monolithic_contact_strategy.cpp
    src/drt_ssi/ssi_monolithic_convcheck_strategies.cpp
    src/drt_ssi/ssi_monolithic_dbc_handler.cpp
    src/drt_ssi/ssi_monolithic_evaluate_OffDiag.cpp
    src/drt_ssi/ssi_monolithic_meshtying_strategy.cpp
    src/drt_ssi/ssi_partitioned.cpp
    src/drt_ssi/ssi_partitioned_1wc.cpp
    src/drt_ssi/ssi_partitioned_2wc.cpp
    src/drt_ssi/ssi_resulttest.cpp
    src/drt_ssi/ssi_utils.cpp
    src/drt_ssi/ssi_str_model_evaluator_partitioned.cpp
    )

set(OBJS_DRT_SSTI
    src/drt_ssti/ssti_algorithm.cpp
    src/drt_ssti/ssti_dyn.cpp
    src/drt_ssti/ssti_monolithic.cpp
    src/drt_ssti/ssti_monolithic_assemble_strategy.cpp
    src/drt_ssti/ssti_monolithic_evaluate_OffDiag.cpp
    src/drt_ssti/ssti_resulttest.cpp
    src/drt_ssti/ssti_utils.cpp
    )

set(OBJS_DRT_IMMERSED_LIB src/drt_immersed_problem/immersed_node.cpp)

set(OBJS_DRT_IMMERSED_PROBLEM
    src/drt_immersed_problem/immersed_base.cpp
    src/drt_immersed_problem/immersed_problem_dyn.cpp
    src/drt_immersed_problem/immersed_partitioned_fsi_dirichletneumann.cpp
    src/drt_immersed_problem/fsi_partitioned_immersed.cpp
    )

set(OBJS_DRT_MOR src/drt_mor/mor_pod.cpp)

set(OBJS_DRT_PASI
    src/drt_pasi/pasi_dyn.cpp
    src/drt_pasi/pasi_partitioned_onewaycoup.cpp
    src/drt_pasi/pasi_partitioned_twowaycoup.cpp
    src/drt_pasi/pasi_partitioned.cpp
    src/drt_pasi/pasi_str_model_evaluator_partitioned.cpp
    src/drt_pasi/pasi_utils.cpp
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
    src/drt_fiber/drt_fiber_utils.cpp
    src/drt_fiber/nodal_fiber_holder.cpp
    )

set(OBJS_DRT_GEOMETRY_PAIR
    src/drt_geometry_pair/geometry_pair_factory.cpp
    src/drt_geometry_pair/geometry_pair.cpp
    src/drt_geometry_pair/geometry_pair_element_faces.cpp
    src/drt_geometry_pair/geometry_pair_evaluation_data_base.cpp
    src/drt_geometry_pair/geometry_pair_line_projection.cpp
    src/drt_geometry_pair/geometry_pair_line_to_3D_evaluation_data.cpp
    src/drt_geometry_pair/geometry_pair_line_to_surface.cpp
    src/drt_geometry_pair/geometry_pair_line_to_surface_gauss_point_projection.cpp
    src/drt_geometry_pair/geometry_pair_line_to_surface_segmentation.cpp
    src/drt_geometry_pair/geometry_pair_line_to_surface_evaluation_data.cpp
    src/drt_geometry_pair/geometry_pair_line_to_volume.cpp
    src/drt_geometry_pair/geometry_pair_line_to_volume_gauss_point_projection.cpp
    src/drt_geometry_pair/geometry_pair_line_to_volume_gauss_point_projection_cross_section.cpp
    src/drt_geometry_pair/geometry_pair_line_to_volume_segmentation.cpp
    src/drt_geometry_pair/geometry_pair_utility_functions.cpp
    )
set(OBJS_DRT_MIXTURE
    src/drt_mixture/growth_evolution_linear_cauchy.H
    src/drt_mixture/growth_evolution.H
    src/drt_mixture/mixture_constituent.cpp
    src/drt_mixture/mixture_constituent_elasthyperbase.cpp
    src/drt_mixture/mixture_constituent_elasthyper.cpp
    src/drt_mixture/mixture_constituent_elasthyper_damage.cpp
    src/drt_mixture/mixture_constituent_elasthyper_elastin_membrane.cpp
    src/drt_mixture/mixture_constituent_muscle_weickenmeier.cpp
    src/drt_mixture/mixture_constituent_remodelfiber_expl.cpp
    src/drt_mixture/mixture_constituent_remodelfiber_impl.cpp
    src/drt_mixture/mixture_constituent_remodelfiber_lib.cpp
    src/drt_mixture/mixture_constituent_remodelfiber_material.cpp
    src/drt_mixture/mixture_constituent_remodelfiber_material_exponential.cpp
    src/drt_mixture/mixture_constituent_remodelfiber_material_exponential_active.cpp
    src/drt_mixture/mixture_growth_strategy.cpp
    src/drt_mixture/mixture_growth_strategy_anisotropic.cpp
    src/drt_mixture/mixture_growth_strategy_isotropic.cpp
    src/drt_mixture/mixture_growth_strategy_stiffness.cpp
    src/drt_mixture/mixture_rule.cpp
    src/drt_mixture/mixture_rule_growthremodel.cpp
    src/drt_mixture/mixture_rule_simple.cpp
    src/drt_mixture/mixture_prestress_strategy.cpp
    src/drt_mixture/mixture_prestress_strategy_constant.cpp
    src/drt_mixture/mixture_prestress_strategy_isocyl.cpp
    src/drt_mixture/mixture_prestress_strategy_iterative.cpp
    src/drt_mixture/remodelfiber.cpp
    )

### LEAGCY
# Define variables that gather the names of all to-be-created libraries. These will be used within
# deprecated_baci_add_library_linked_with_everything() to link every target with everything.

set(ELEMENT_LIBRARIES
    ele_ale
    ele_beam3
    ele_rigidsphere
    ele_bele
    ele_discsh3
    ele_fluid
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
    ele_elemag
    )

set(MULTIFIELD_LIBRARIES
    alg_elch
    alg_fsi
    alg_fbi
    alg_fsi_xfem
    alg_loma
    alg_two_phase
    alg_tsi
    alg_sti
    alg_wear
    alg_poroelast
    alg_poromultiphase
    alg_immersed_problem
    alg_fpsi
    alg_fs3i
    alg_ssi
    alg_ssti
    alg_ehl
    alg_pasi
    )

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
    alg_statinvanalysis
    alg_thermo
    alg_contact
    alg_contact_aug
    alg_particle_algorithm
    alg_particle_engine
    alg_particle_interaction
    alg_particle_rigidbody
    alg_particle_wall
    alg_elemag
    alg_tutorial
    )

set(SUPPORT_LIBRARIES
    drt_lib_xfem_manag
    drt_lib
    drt_lib_artnet
    drt_lib_browniandyn
    drt_lib_redairways
    drt_lib_mortar
    drt_lib_constraint
    drt_lib_contact
    drt_lib_beamcontact
    drt_lib_beaminteraction
    drt_levelset
    drt_lib_poromulti_art_coupl
    # Libraries needed by elements
    drt_lib_surfstress
    drt_lib_xfem
    drt_lib_cut
    drt_lib_inpar
    drt_lib_io
    drt_io
    drt_lib_contact_constitutivelaw
    drt_lib_nurbs
    drt_lib_binstrategy
    drt_lib_particleengine
    drt_lib_statinv_ana
    drt_lib_immersed
    drt_lib_fiber
    drt_lib_geometric_search
    drt_lib_geometry_pair
    )

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
    solver
    linalg
    drt_comm
    headers
    drt_timestepping
    )

###

###------------------------------------------------------------------ Binaries

# hide libraries

if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/src)
endif(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

add_library(headers INTERFACE)
target_include_directories(headers INTERFACE ${PROJECT_SOURCE_DIR}/src/headers)

add_library(drt_timestepping INTERFACE)
target_include_directories(drt_timestepping INTERFACE ${PROJECT_SOURCE_DIR}/src/drt_timestepping)

# Compile element libraries first. There are files that take time to compile.
deprecated_baci_add_library_linked_with_everything(ele_rigidsphere ${OBJS_DRT_RIGIDSPHERE_LIB})
deprecated_baci_add_library_linked_with_everything(ele_s8 ${OBJS_DRT_S8_LIB})
deprecated_baci_add_library_linked_with_everything(ele_membrane ${OBJS_DRT_MEMBRANE_LIB})
deprecated_baci_add_library_linked_with_everything(
  ele_porofluidmultiphase ${OBJS_DRT_POROFLUIDMULTIPHASE_LIB}
  )
deprecated_baci_add_library_linked_with_everything(ele_scatra ${OBJS_DRT_SCATRA_LIB})
deprecated_baci_add_library_linked_with_everything(ele_sp3 ${OBJS_DRT_SP3_LIB})
deprecated_baci_add_library_linked_with_everything(ele_so ${OBJS_SOH8_LIB})
deprecated_baci_add_library_linked_with_everything(ele_t3 ${OBJS_DRT_T3_LIB})
deprecated_baci_add_library_linked_with_everything(ele_t3cl ${OBJS_DRT_T3CL_LIB})
deprecated_baci_add_library_linked_with_everything(ele_to3 ${OBJS_DRT_TO3_LIB})
deprecated_baci_add_library_linked_with_everything(ele_w1 ${OBJS_DRT_W1_LIB})

# algorithm libraries that depend on adapters
deprecated_baci_add_library_linked_with_everything(alg_fbi ${OBJS_DRT_FBI})
deprecated_baci_add_library_linked_with_everything(alg_fsi ${OBJS_DRT_FSI})
deprecated_baci_add_library_linked_with_everything(
  alg_immersed_problem ${OBJS_DRT_IMMERSED_PROBLEM}
  )
deprecated_baci_add_library_linked_with_everything(alg_fsi_xfem ${OBJS_DRT_FSI_XFEM})
deprecated_baci_add_library_linked_with_everything(alg_two_phase ${OBJS_DRT_TWO_PHASE})
deprecated_baci_add_library_linked_with_everything(alg_sti ${OBJS_DRT_STI})
deprecated_baci_add_library_linked_with_everything(alg_wear ${OBJS_DRT_WEAR})
deprecated_baci_add_library_linked_with_everything(alg_poroelast ${OBJS_DRT_POROELAST})
deprecated_baci_add_library_linked_with_everything(alg_poromultiphase ${OBJS_DRT_POROMULTIPHASE})
deprecated_baci_add_library_linked_with_everything(alg_fpsi ${OBJS_DRT_FPSI})
deprecated_baci_add_library_linked_with_everything(alg_fs3i ${OBJS_DRT_FS3I})
deprecated_baci_add_library_linked_with_everything(alg_ssi ${OBJS_DRT_SSI})
deprecated_baci_add_library_linked_with_everything(alg_ssti ${OBJS_DRT_SSTI})
deprecated_baci_add_library_linked_with_everything(alg_pasi ${OBJS_DRT_PASI})

# algorithm libraries the adapters depend uo
deprecated_baci_add_library_linked_with_everything(alg_statinvanalysis ${OBJS_STATINV_ANALYSIS})
deprecated_baci_add_library_linked_with_everything(alg_redairways ${OBJS_DRT_RED_AIRWAYS})
deprecated_baci_add_library_linked_with_everything(
  alg_porofluidmultiphase ${OBJS_DRT_POROFLUIDMULTIPHASE}
  )
deprecated_baci_add_library_linked_with_everything(alg_scatra ${OBJS_DRT_SCATRA})
deprecated_baci_add_library_linked_with_everything(alg_structure ${OBJS_DRT_STRUCTURE})
deprecated_baci_add_library_linked_with_everything(alg_xstructure ${OBJS_DRT_STRUCTURE_XSTRUCTURE})
deprecated_baci_add_library_linked_with_everything(
  alg_particle_algorithm ${OBJS_DRT_PARTICLE_ALGORITHM}
  )
deprecated_baci_add_library_linked_with_everything(alg_particle_engine ${OBJS_DRT_PARTICLE_ENGINE})
deprecated_baci_add_library_linked_with_everything(
  alg_particle_interaction ${OBJS_DRT_PARTICLE_INTERACTION}
  )
deprecated_baci_add_library_linked_with_everything(
  alg_particle_rigidbody ${OBJS_DRT_PARTICLE_RIGIDBODY}
  )
deprecated_baci_add_library_linked_with_everything(alg_particle_wall ${OBJS_DRT_PARTICLE_WALL})
deprecated_baci_add_library_linked_with_everything(alg_tutorial ${OBJS_DRT_TUTORIAL})

# general framework
deprecated_baci_add_library_linked_with_everything(drt_nurbs ${OBJS_DRT_NURBS_DISCRET})
deprecated_baci_add_library_linked_with_everything(drt_mortar ${OBJS_DRT_MORTAR})
deprecated_baci_add_library_linked_with_everything(drt_volmortar ${OBJS_DRT_VOLMORTAR})
deprecated_baci_add_library_linked_with_everything(drt_mor ${OBJS_DRT_MOR})

# material definitions
deprecated_baci_add_library_linked_with_everything(drt_mat ${OBJS_DRT_MAT_LIB} ${OBJS_DRT_MIXTURE})
deprecated_baci_add_library_linked_with_everything(
  drt_micromat ${OBJS_DRT_MAT} ${OBJS_DRT_STRU_MULTI}
  )
deprecated_baci_add_library_linked_with_everything(drt_scatra_multi ${OBJS_DRT_MAT_SCATRA_MULTI})

# special library that demands all parobject types
deprecated_baci_add_library_linked_with_everything(drt_lib_register ${OBJS_DRT_REGISTER})

# libraries that are used by baci and by filters
deprecated_baci_add_library_linked_with_everything(drt_lib_inpar ${OBJS_DRT_INPAR})
deprecated_baci_add_library_linked_with_everything(
  drt_lib_poromulti_art_coupl ${OBJS_DRT_POROMULTIPHASE_ART_COUPL}
  )
deprecated_baci_add_library_linked_with_everything(drt_lib_xfem ${OBJS_DRT_XFEM})
deprecated_baci_add_library_linked_with_everything(drt_lib_mortar ${OBJS_DRT_MORTAR_LIB})
deprecated_baci_add_library_linked_with_everything(drt_lib_nurbs ${OBJS_DRT_NURBS_DISCRET_LIB})

deprecated_baci_add_library_linked_with_everything(
  drt_lib_particleengine ${OBJS_DRT_PARTICLE_ENGINE_LIB}
  )
deprecated_baci_add_library_linked_with_everything(drt_lib_surfstress ${OBJS_DRT_SURFSTRESS})
deprecated_baci_add_library_linked_with_everything(drt_lib_utils ${OBJS_DRT_UTILS})
deprecated_baci_add_library_linked_with_everything(drt_lib_redairways ${OBJS_DRT_RED_AIRWAYS_LIB})
deprecated_baci_add_library_linked_with_everything(drt_lib_statinv_ana ${OBJS_STATINV_ANALYSIS_LIB})
deprecated_baci_add_library_linked_with_everything(drt_lib_immersed ${OBJS_DRT_IMMERSED_LIB})
deprecated_baci_add_library_linked_with_everything(drt_lib_xfem_manag ${OBJS_DRT_XFEM_MANAG_LIB})
deprecated_baci_add_library_linked_with_everything(drt_lib_fiber ${OBJS_DRT_FIBER_LIB})
deprecated_baci_add_library_linked_with_everything(drt_lib_geometry_pair ${OBJS_DRT_GEOMETRY_PAIR})
deprecated_baci_add_library_linked_with_everything(alg_mor ${OBJS_DRT_MOR})

# filter framework
deprecated_baci_add_library_linked_with_everything(post_drt_common ${OBJS_POST_DRT_COMMON})
deprecated_baci_add_library_linked_with_everything(post_drt_ensight_lib ${OBJS_POST_DRT_ENSIGHT})
deprecated_baci_add_library_linked_with_everything(post_drt_vtk_lib ${OBJS_POST_DRT_VTK})

target_link_libraries(post_drt_vtk_lib PUBLIC post_drt_common)
target_link_libraries(post_drt_ensight_lib PUBLIC post_drt_common)

# continuation method
deprecated_baci_add_library_linked_with_everything(locacontinuation ${OBJS_LOCA_CONTINUATION})

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
