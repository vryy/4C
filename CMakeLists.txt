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
    )
set(OBJS_POST_DRT_MONITOR src/post_monitor/post_drt_monitor.cpp)
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

set(OBJS_DRT_REGISTER src/drt_lib/drt_parobjectregister.cpp)

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

set(OBJS_DRT_MOR src/drt_mor/mor_pod.cpp)

set(OBJS_DRT_PASI
    src/drt_pasi/pasi_dyn.cpp
    src/drt_pasi/pasi_partitioned_onewaycoup.cpp
    src/drt_pasi/pasi_partitioned_twowaycoup.cpp
    src/drt_pasi/pasi_partitioned.cpp
    src/drt_pasi/pasi_str_model_evaluator_partitioned.cpp
    src/drt_pasi/pasi_utils.cpp
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
    alg_poromultiphase_scatra
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
    alg_lubrication
    alg_porofluidmultiphase
    alg_scatra
    alg_structure
    alg_structure_new
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
    drt_red_airways
    drt_lib_mortar
    drt_lib_constraint
    drt_lib_contact
    drt_lib_beamcontact
    drt_lib_beaminteraction
    drt_levelset
    drt_lib_poromulti_art_coupl
    # Libraries needed by elements
    drt_lib_struct_ale
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
    drt_stru_multi
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

# Compile element libraries first. There are files that take time to compile.
deprecated_baci_add_library_linked_with_everything(ele_rigidsphere ${OBJS_DRT_RIGIDSPHERE_LIB})
deprecated_baci_add_library_linked_with_everything(ele_s8 ${OBJS_DRT_S8_LIB})
deprecated_baci_add_library_linked_with_everything(ele_membrane ${OBJS_DRT_MEMBRANE_LIB})

deprecated_baci_add_library_linked_with_everything(ele_scatra ${OBJS_DRT_SCATRA_LIB})

# algorithm libraries that depend on adapters
deprecated_baci_add_library_linked_with_everything(alg_ssi ${OBJS_DRT_SSI})
deprecated_baci_add_library_linked_with_everything(alg_sti ${OBJS_DRT_STI})
deprecated_baci_add_library_linked_with_everything(alg_ssti ${OBJS_DRT_SSTI})
deprecated_baci_add_library_linked_with_everything(alg_pasi ${OBJS_DRT_PASI})

# algorithm libraries the adapters depend uo
deprecated_baci_add_library_linked_with_everything(alg_statinvanalysis ${OBJS_STATINV_ANALYSIS})
deprecated_baci_add_library_linked_with_everything(alg_scatra ${OBJS_DRT_SCATRA})
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

# general framework
deprecated_baci_add_library_linked_with_everything(drt_nurbs ${OBJS_DRT_NURBS_DISCRET})
deprecated_baci_add_library_linked_with_everything(drt_mor ${OBJS_DRT_MOR})

# material definitions
deprecated_baci_add_library_linked_with_everything(drt_mat ${OBJS_DRT_MAT_LIB} ${OBJS_DRT_MIXTURE})
deprecated_baci_add_library_linked_with_everything(drt_micromat ${OBJS_DRT_MAT})
deprecated_baci_add_library_linked_with_everything(drt_scatra_multi ${OBJS_DRT_MAT_SCATRA_MULTI})

# special library that demands all parobject types
deprecated_baci_add_library_linked_with_everything(drt_lib_register ${OBJS_DRT_REGISTER})

# libraries that are used by baci and by filters
deprecated_baci_add_library_linked_with_everything(drt_lib_inpar ${OBJS_DRT_INPAR})
deprecated_baci_add_library_linked_with_everything(drt_lib_nurbs ${OBJS_DRT_NURBS_DISCRET_LIB})

deprecated_baci_add_library_linked_with_everything(
  drt_lib_particleengine ${OBJS_DRT_PARTICLE_ENGINE_LIB}
  )
deprecated_baci_add_library_linked_with_everything(drt_lib_utils ${OBJS_DRT_UTILS})
deprecated_baci_add_library_linked_with_everything(drt_lib_statinv_ana ${OBJS_STATINV_ANALYSIS_LIB})
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
