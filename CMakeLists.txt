cmake_minimum_required(VERSION 2.8)
#project(baci CXX C Fortran)
project(baci CXX C Fortran)

# The version number.
set (baci_VERSION_MAJOR 1)
set (baci_VERSION_MINOR 0)

# configure a header file to pass some of the CMake settings
# to the source code
#configure_file (
#  "${PROJECT_SOURCE_DIR}/baci_config.h.in"
#  "${PROJECT_BINARY_DIR}/baci_config.h"
#  )
configure_file (
  "${PROJECT_SOURCE_DIR}/src/headers/compile_settings.h.in"
  "${PROJECT_BINARY_DIR}/src/headers/compile_settings.h"
  )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}/src/headers")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules")

set(LIBRARIES "")

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)
find_package(Trilinos REQUIRED)

#find_package(Blitz)

SET(Boost_USE_STATIC_LIBS ON) 
SET(Boost_USE_MULTITHREAD OFF) 
#FIND_PACKAGE( Boost COMPONENTS date_time filesystem iostreams program_options )
FIND_PACKAGE( Boost COMPONENTS program_options )

if (MPI_FOUND)
    #include_directories(${MPI_INCLUDE_PATH})
    SET(CMAKE_C_COMPILER mpicc)
    SET(CMAKE_CXX_COMPILER mpic++)
endif(MPI_FOUND)

# trilinos libraries
if (Trilinos_FOUND)
  include_directories(${Trilinos_INCLUDE_DIRS})
  link_directories(${Trilinos_LIBRARY_DIRS})
  set(LIBRARIES
    ${LIBRARIES}
    #${Trilinos_LIBRARIES}
    trilinoscouplings
    noxepetra
    nox
    ifpack
    aztecoo
    amesos
    ml
    zoltan
    epetraext
    epetra
    teuchos
    #stk_mesh_fem stk_mesh_base stk_util_use_cases stk_util_parallel stk_util_diag stk_util_env stk_util_util
    )
endif (Trilinos_FOUND)

find_package(Umfpack)
if (Umfpack_FOUND)
  #message("Umfpack found")
  include_directories(${Umfpack_INCLUDE_DIRS})
  link_directories(${Umfpack_LIBRARY_DIRS})
  set(LIBRARIES
    ${LIBRARIES}
    ${Umfpack_LIBRARIES}
    )
else (Umfpack_FOUND)
  #message("Umfpack not found")
  add_definitions("-DHAVENOT_UMFPACK")
endif (Umfpack_FOUND)

find_package(SuperLU)
if (SuperLU_FOUND)
  include_directories(${SuperLU_INCLUDE_DIRS})
  link_directories(${SuperLU_LIBRARY_DIRS})
  set(LIBRARIES
    ${LIBRARIES}
    ${SuperLU_LIBRARIES}
    )
else (SuperLU_FOUND)
  add_definitions("-DHAVENOT_SUPERLU")
endif (SuperLU_FOUND)

find_package(PARMETIS)
if (PARMETIS_FOUND)
  include_directories(${PARMETIS_INCLUDE_DIRS})
  link_directories(${PARMETIS_LIBRARY_DIRS})
  set(LIBRARIES
    ${LIBRARIES}
    ${PARMETIS_LIBRARIES}
    )
endif (PARMETIS_FOUND)

if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  link_directories(${Boost_LIBRARY_DIRS})
  set(LIBRARIES
    ${LIBRARIES}
    ${Boost_LIBRARIES}
    )
endif(Boost_FOUND)

if (LAPACK_FOUND)
  include_directories(${LAPACK_INCLUDE_DIRS})
  link_directories(${LAPACK_LIBRARY_DIRS})
  set(LIBRARIES
    ${LIBRARIES}
    ${LAPACK_LIBRARIES}
    )
endif(LAPACK_FOUND)

if (BLAS_FOUND)
  include_directories(${BLAS_INCLUDE_DIRS})
  link_directories(${BLAS_LIBRARY_DIRS})
  set(LIBRARIES
    ${LIBRARIES}
    ${BLAS_LIBRARIES}
    )
endif(BLAS_FOUND)

add_definitions("-DLINUX_MUENCH")
add_definitions("-DPARALLEL")
add_definitions("-DLOOPL=15")
add_definitions("-DCCADISCRET")
add_definitions("-DTRILINOS_PACKAGE")
add_definitions("-DD_FSI")

### compile time options

option (DEBUG "Turn on special debug code" ON)
if (DEBUG)
  add_definitions("-DDEBUG")
endif (DEBUG)

set( adapter_source
src/drt_adapter/adapter_algorithmbase.cpp
src/drt_adapter/adapter_coupling.cpp
src/drt_adapter/adapter_coupling_mortar.cpp
src/drt_adapter/adapter_fluid_ale.cpp
src/drt_adapter/adapter_fluid_base_algorithm.cpp
src/drt_adapter/adapter_fluid_combust.cpp
src/drt_adapter/adapter_fluid.cpp
src/drt_adapter/adapter_fluid_genalpha.cpp
src/drt_adapter/adapter_fluid_impl.cpp
src/drt_adapter/adapter_fluid_lung.cpp
src/drt_adapter/adapter_fluid_moving_boundary.cpp
src/drt_adapter/adapter_fluid_projection.cpp
src/drt_adapter/adapter_fluid_xfem.cpp
src/drt_adapter/adapter_scatra_base_algorithm.cpp
src/drt_adapter/adapter_scatra_fluid_ale_coupling_algo.cpp
src/drt_adapter/adapter_scatra_fluid_coupling_algorithm.cpp
src/drt_adapter/adapter_structure_cmtstrugenalpha.cpp
src/drt_adapter/adapter_structure_constr_merged.cpp
src/drt_adapter/adapter_structure.cpp
src/drt_adapter/adapter_structure_lung.cpp
src/drt_adapter/adapter_structure_strugenalpha.cpp
src/drt_adapter/adapter_structure_timint.cpp
src/drt_adapter/adapter_structure_wrapper.cpp
src/drt_adapter/adapter_thermo.cpp
src/drt_adapter/adapter_thermo_timint.cpp
src/drt_adapter/adapter_xfluid_impl.cpp
)
set( ale_source
src/drt_ale/ale.cpp
src/drt_ale/ale_dyn.cpp
src/drt_ale/ale_laplace.cpp
src/drt_ale/ale_lin.cpp
src/drt_ale/ale_resulttest.cpp
src/drt_ale/ale_springs.cpp
src/drt_ale/ale_springs_fixed_ref.cpp
src/drt_ale/ale_utils_mapextractor.cpp
)
set( ale2_source
src/drt_ale2/ale2.cpp
src/drt_ale2/ale2_evaluate.cpp
src/drt_ale2/ale2_input.cpp
src/drt_ale2/ale2_line.cpp
src/drt_ale2/ale2_line_evaluate.cpp
src/drt_ale2/ale2_nurbs.cpp
)
set( ale3_source
src/drt_ale3/ale3.cpp
src/drt_ale3/ale3_evaluate.cpp
src/drt_ale3/ale3_input.cpp
src/drt_ale3/ale3_nurbs.cpp
src/drt_ale3/ale3_surface.cpp
src/drt_ale3/ale3_surface_evaluate.cpp
)
set( bele3_source
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
set( combust_element_source
src/drt_combust/combust3.cpp
src/drt_combust/combust3_evaluate.cpp
src/drt_combust/combust3_input.cpp
src/drt_combust/combust3_line.cpp
src/drt_combust/combust3_line_evaluate.cpp
src/drt_combust/combust3_surface.cpp
src/drt_combust/combust3_surface_evaluate.cpp
src/drt_combust/combust3_sysmat.cpp
)
set( combust_source
src/drt_combust/combust_algorithm.cpp
src/drt_combust/combust_dyn.cpp
src/drt_combust/combust_flamefront.cpp
src/drt_combust/combust_fluidimplicitintegration.cpp
src/drt_combust/combust_fluidresulttest.cpp
src/drt_combust/combust_interface.cpp
src/drt_combust/combust_refinementcell.cpp
src/drt_combust/combust_reinitializer.cpp
src/drt_combust/combust_utils.cpp
)
set( constraint_element_source
src/drt_constraint/constraint_element2.cpp
src/drt_constraint/constraint_element2_evaluate.cpp
src/drt_constraint/constraint_element3.cpp
src/drt_constraint/constraint_element3_evaluate.cpp
)
set( constraint_source
src/drt_constraint/constraint.cpp
src/drt_constraint/constraintdofset.cpp
src/drt_constraint/constraint_manager.cpp
src/drt_constraint/constraintsolver.cpp
src/drt_constraint/monitor.cpp
src/drt_constraint/mpcdofset.cpp
src/drt_constraint/multipointconstraint2.cpp
src/drt_constraint/multipointconstraint3.cpp
src/drt_constraint/multipointconstraint.cpp
)
set( contact_source
src/drt_contact/beam3contact.cpp
src/drt_contact/beam3contact_manager.cpp
src/drt_contact/beam3contactstrugenalpha.cpp
src/drt_contact/contactstrugenalpha.cpp
src/drt_contact/contactstrugenalpha_lsnewton.cpp
src/drt_contact/drt_cdofset.cpp
src/drt_contact/drt_celement.cpp
src/drt_contact/drt_celement_shapefct.cpp
src/drt_contact/drt_cnode.cpp
src/drt_contact/drt_cnodedatacontainer.cpp
src/drt_contact/drt_contact_abstract_strategy.cpp
src/drt_contact/drt_contact_binarytree.cpp
src/drt_contact/drt_contact_binarytree_self.cpp
src/drt_contact/drt_contact_coupling2d.cpp
src/drt_contact/drt_contact_coupling3d_classes.cpp
src/drt_contact/drt_contact_coupling3d.cpp
src/drt_contact/drt_contact_coupling3d_tools.cpp
src/drt_contact/drt_contact_integrator.cpp
src/drt_contact/drt_contact_integrator_tools.cpp
src/drt_contact/drt_contact_interface.cpp
src/drt_contact/drt_contact_interface_tools.cpp
src/drt_contact/drt_contact_lagrange_strategy.cpp
src/drt_contact/drt_contact_manager_base.cpp
src/drt_contact/drt_contact_manager.cpp
src/drt_contact/drt_contact_penalty_strategy.cpp
src/drt_contact/drt_contact_projector.cpp
)
set( contactnew_source
src/drt_contactnew/contact_abstract_strategy.cpp
src/drt_contactnew/contact_coupling2d.cpp
src/drt_contactnew/contact_coupling3d.cpp
src/drt_contactnew/contact_element.cpp
src/drt_contactnew/contact_integrator.cpp
src/drt_contactnew/contact_interface.cpp
src/drt_contactnew/contact_lagrange_strategy.cpp
src/drt_contactnew/contact_manager.cpp
src/drt_contactnew/contact_node.cpp
src/drt_contactnew/contact_penalty_strategy.cpp
src/drt_contactnew/friction_node.cpp
src/drt_contactnew/meshtying_abstract_strategy.cpp
src/drt_contactnew/meshtying_lagrange_strategy.cpp
src/drt_contactnew/meshtying_manager.cpp
src/drt_contactnew/meshtying_penalty_strategy.cpp
src/drt_contactnew/selfcontact_binarytree.cpp
src/drt_contactnew/strugenalpha_cmt.cpp
)
set( elch_source
src/drt_elch/elch_algorithm.cpp
src/drt_elch/elch_dyn.cpp
src/drt_elch/elch_moving_boundary_algorithm.cpp
)
set( fem_general_source
src/drt_fem_general/drt_utils_boundary_integration.cpp
src/drt_fem_general/drt_utils_bspline.cpp
src/drt_fem_general/drt_utils_dgfem_basisfunctions.cpp
src/drt_fem_general/drt_utils_fem_shapefunctions.cpp
src/drt_fem_general/drt_utils_integration.cpp
src/drt_fem_general/drt_utils_local_connectivity_matrices.cpp
src/drt_fem_general/drt_utils_nurbs_shapefunctions.cpp
)
set( fluid_element_source
src/drt_fluid/time_integration_element.cpp
)
set( fluid_source
src/drt_fluid/drt_pbcdofset.cpp
src/drt_fluid/drt_periodicbc.cpp
src/drt_fluid/dyn_smag.cpp
src/drt_fluid/fluid_coupling_red_models.cpp
src/drt_fluid/fluid_dyn_nln_drt.cpp
src/drt_fluid/fluid_genalpha_integration.cpp
src/drt_fluid/fluidimpedancecondition.cpp
src/drt_fluid/fluidimplicitintegration.cpp
src/drt_fluid/fluid_projectionmethod.cpp
src/drt_fluid/fluidresulttest.cpp
src/drt_fluid/fluid_rotsym_periodicbc_utils.cpp
src/drt_fluid/fluid_utils.cpp
src/drt_fluid/fluid_utils_mapextractor.cpp
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
)
set( fsi_source
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
src/drt_fsi/fsi_monolithiclagrange.cpp
src/drt_fsi/fsi_monolithic_linearsystem.cpp
src/drt_fsi/fsi_monolithicoverlap.cpp
src/drt_fsi/fsi_monolithicstructuresplit.cpp
src/drt_fsi/fsi_monolithicxfem.cpp
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
src/drt_fsi/fsi_partitioned.cpp
src/drt_fsi/fsi_partitionedmonolithic.cpp
src/drt_fsi/fsi_robin.cpp
src/drt_fsi/fsi_robinneumann.cpp
src/drt_fsi/fsi_statustest.cpp
src/drt_fsi/fsi_structureale.cpp
src/drt_fsi/fsi_utils.cpp
src/drt_fsi/fs_monolithic.cpp
)
set( geometry_source
src/drt_geometry/integrationcell.cpp
src/drt_geometry/intersection.cpp
src/drt_geometry/intersection_interfacepoint.cpp
src/drt_geometry/intersection_math.cpp
src/drt_geometry/intersection_service.cpp
src/drt_geometry/position_array.cpp
src/drt_geometry/searchtree.cpp
src/drt_geometry/searchtree_geometry_service.cpp
src/drt_geometry/searchtree_nearestobject.cpp
)
set( inpar_source
src/drt_inpar/drt_boolifyparameters.cpp
src/drt_inpar/drt_validconditions.cpp
src/drt_inpar/drt_validmaterials.cpp
src/drt_inpar/drt_validparameters.cpp
)
set( inv_analysis_source
src/drt_inv_analysis/gen_inv_analysis.cpp
src/drt_inv_analysis/inv_analysis.cpp
)
set( io_source
src/drt_io/io_control.cpp
src/drt_io/io.cpp
src/drt_io/io_gmsh.cpp
src/drt_io/io_gmsh_xfem_extension.cpp
src/drt_io/io_hdf.cpp
)
set( lib_element_source
src/drt_lib/drt_utils.cpp
)
set( lib_source
src/drt_lib/distribute_grids.cpp
src/drt_lib/drt_condition.cpp
src/drt_lib/drt_conditiondefinition.cpp
src/drt_lib/drt_condition_selector.cpp
src/drt_lib/drt_condition_utils.cpp
src/drt_lib/drt_container.cpp
src/drt_lib/drt_discret_conditions.cpp
src/drt_lib/drt_discret.cpp
src/drt_lib/drt_discret_evaluate.cpp
src/drt_lib/drt_discret_fillcomplete.cpp
src/drt_lib/drt_discret_partition.cpp
src/drt_lib/drt_discret_utils.cpp
src/drt_lib/drt_dofset_base.cpp
src/drt_lib/drt_dofset.cpp
src/drt_lib/drt_dofset_proxy.cpp
src/drt_lib/drt_dserror.cpp
src/drt_lib/drt_element.cpp
src/drt_lib/drt_elementdefinition.cpp
src/drt_lib/drt_elementreader.cpp
src/drt_lib/drt_elementregister.cpp
src/drt_lib/drt_exporter.cpp
src/drt_lib/drt_function.cpp
src/drt_lib/drt_globalproblem.cpp
src/drt_lib/drt_init_control.cpp
src/drt_lib/drt_inputreader.cpp
src/drt_lib/drt_linedefinition.cpp
src/drt_lib/drt_locsys.cpp
src/drt_lib/drt_materialdefinition.cpp
src/drt_lib/drt_node.cpp
src/drt_lib/drt_nodematchingoctree.cpp
src/drt_lib/drt_nodereader.cpp
src/drt_lib/drt_parobject.cpp
src/drt_lib/drt_parser.cpp
src/drt_lib/drt_resulttest.cpp
src/drt_lib/drt_timecurve.cpp
src/drt_lib/drt_utils_createdis.cpp
src/drt_lib/drt_utils_metis.cpp
src/drt_lib/drt_utils_parmetis.cpp
src/drt_lib/global_inp_control2.cpp
src/drt_lib/linalg_ana.cpp
src/drt_lib/linalg_blocksparsematrix.cpp
src/drt_lib/linalg_condest.cpp
src/drt_lib/linalg_downwindmatrix.cpp
src/drt_lib/linalg_krylov_projector.cpp
src/drt_lib/linalg_mapextractor.cpp
src/drt_lib/linalg_mlapi_operator.cpp
src/drt_lib/linalg_precond.cpp
src/drt_lib/linalg_projected_operator.cpp
src/drt_lib/linalg_projected_precond.cpp
src/drt_lib/linalg_serialdensematrix.cpp
src/drt_lib/linalg_serialdensevector.cpp
src/drt_lib/linalg_solver.cpp
src/drt_lib/linalg_solver_spooles.cpp
src/drt_lib/linalg_sparsematrix.cpp
src/drt_lib/linalg_utils.cpp
src/drt_lib/simpler_operator.cpp
)
set( loma_source
src/drt_loma/loma_algorithm.cpp
src/drt_loma/loma_dyn.cpp
)
set( mat_source
src/drt_mat/material.cpp
src/drt_mat/matpar_bundle.cpp
src/drt_mat/matpar_material.cpp
src/drt_mat/matpar_parameter.cpp
src/drt_mat/aaaneohooke.cpp
src/drt_mat/aaaraghavanvorp_damage.cpp
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
src/drt_mat/contchainnetw.cpp
src/drt_mat/elasthyper.cpp
src/drt_mat/ferech_pv.cpp
src/drt_mat/fourieriso.cpp
src/drt_mat/holzapfelcardiovascular.cpp
src/drt_mat/humphreycardiovascular.cpp
src/drt_mat/ion.cpp
src/drt_mat/itskov.cpp
src/drt_mat/logneohooke.cpp
src/drt_mat/lung_ogden.cpp
src/drt_mat/lung_penalty.cpp
src/drt_mat/matlist.cpp
src/drt_mat/micromaterial.cpp
src/drt_mat/micromaterial_evaluate.cpp
src/drt_mat/micromaterialgp_static.cpp
src/drt_mat/mixfrac.cpp
src/drt_mat/modpowerlaw.cpp
src/drt_mat/mooneyrivlin.cpp
src/drt_mat/neohooke.cpp
src/drt_mat/newtonianfluid.cpp
src/drt_mat/plasticneohooke.cpp
src/drt_mat/protein.cpp
src/drt_mat/scatra_mat.cpp
src/drt_mat/stvenantkirchhoff.cpp
src/drt_mat/sutherland.cpp
src/drt_mat/viscoanisotropic.cpp
src/drt_mat/visconeohooke.cpp
src/drt_mat/yeoh.cpp
)
set( matelast_source
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
src/drt_matelast/elast_vologden.cpp
src/drt_matelast/elast_volpenalty.cpp
src/drt_matelast/elast_volsussmanbathe.cpp
)
set( mortar_source
src/drt_mortar/mortar_binarytree.cpp
src/drt_mortar/mortar_coupling2d.cpp
src/drt_mortar/mortar_coupling3d_classes.cpp
src/drt_mortar/mortar_coupling3d.cpp
src/drt_mortar/mortar_defines.cpp
src/drt_mortar/mortar_dofset.cpp
src/drt_mortar/mortar_element.cpp
src/drt_mortar/mortar_element_shapefct.cpp
src/drt_mortar/mortar_integrator.cpp
src/drt_mortar/mortar_interface.cpp
src/drt_mortar/mortar_manager_base.cpp
src/drt_mortar/mortar_node.cpp
src/drt_mortar/mortar_projector.cpp
src/drt_mortar/mortar_strategy_base.cpp
)
set( nurbs_discret_source
src/drt_nurbs_discret/drt_apply_nurbs_initial_condition.cpp
src/drt_nurbs_discret/drt_control_point.cpp
src/drt_nurbs_discret/drt_knotvector.cpp
src/drt_nurbs_discret/drt_nurbs_discret.cpp
)
set( potential_source
src/drt_potential/drt_potential_container.cpp
src/drt_potential/drt_potential.cpp
src/drt_potential/drt_potential_dofset.cpp
src/drt_potential/drt_potential_manager.cpp
src/drt_potential/drt_potential_surface.cpp
src/drt_potential/drt_potential_volume.cpp
)
set( scatra_source
src/drt_scatra/passive_scatra_algorithm.cpp
src/drt_scatra/scatra_dyn.cpp
src/drt_scatra/scatra_ele_boundary_impl.cpp
src/drt_scatra/scatra_ele_impl.cpp
src/drt_scatra/scatra_ele_impl_utils.cpp
src/drt_scatra/scatra_element_boundary_evaluate.cpp
src/drt_scatra/scatra_element.cpp
src/drt_scatra/scatra_element_evaluate.cpp
src/drt_scatra/scatra_element_input.cpp
src/drt_scatra/scatra_resulttest.cpp
src/drt_scatra/scatra_timint_bdf2.cpp
src/drt_scatra/scatra_timint_genalpha.cpp
src/drt_scatra/scatra_timint_implicit.cpp
src/drt_scatra/scatra_timint_implicit_service.cpp
src/drt_scatra/scatra_timint_ost.cpp
src/drt_scatra/scatra_timint_stat.cpp
src/drt_scatra/scatra_utils.cpp
)
set( statmech_source
src/drt_statmech/bromotion_manager.cpp
src/drt_statmech/statmech_manager.cpp
src/drt_statmech/statmech_time.cpp
)
set( structure_source
src/drt_structure/str_aux.cpp
src/drt_structure/str_invanalysis.cpp
src/drt_structure/strtimada.cpp
src/drt_structure/strtimada_create.cpp
src/drt_structure/strtimada_zienxie.cpp
src/drt_structure/strtimint_ab2.cpp
src/drt_structure/strtimint.cpp
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
src/drt_structure/structure_utils_mapextractor.cpp
src/drt_structure/strudyn_direct.cpp
src/drt_structure/stru_dyn_nln_drt.cpp
src/drt_structure/strugenalpha.cpp
src/drt_structure/strugenalpha_lsnewton.cpp
src/drt_structure/stru_resulttest.cpp
src/drt_structure/stru_static_drt.cpp
src/structure/stru_static_lin.c
)
set( stru_multi_source
src/drt_stru_multi/convert_mat.cpp
src/drt_stru_multi/microstatic.cpp
)
set( surfstress_source
src/drt_surfstress/drt_surfstress_manager.cpp
)
set( thermo_main_source
src/drt_thermo/thr_dyn.cpp
)
set( thermo_timint_source
src/drt_thermo/thr_resulttest.cpp
src/drt_thermo/thr_aux.cpp
src/drt_thermo/thrtimint.cpp
src/drt_thermo/thrtimint_genalpha.cpp
src/drt_thermo/thrtimint_impl.cpp
src/drt_thermo/thrtimint_ost.cpp
src/drt_thermo/thrtimint_statics.cpp
)
set( tsi_source
src/drt_tsi/tsi_algorithm.cpp
src/drt_tsi/tsi_dyn.cpp
src/drt_tsi/tsi_utils.cpp
)
set( xdiff3_source
src/drt_xdiff3/xdiff3.cpp
src/drt_xdiff3/xdiff3_evaluate.cpp
src/drt_xdiff3/xdiff3_input.cpp
src/drt_xdiff3/xdiff3_line.cpp
src/drt_xdiff3/xdiff3_line_evaluate.cpp
src/drt_xdiff3/xdiff3_surface.cpp
src/drt_xdiff3/xdiff3_surface_evaluate.cpp
src/drt_xdiff3/xdiff3_sysmat4.cpp
)
set( xfem_source
src/drt_xfem/dof_distribution_switcher.cpp
src/drt_xfem/dof_management.cpp
src/drt_xfem/dof_management_element.cpp
src/drt_xfem/enrichment.cpp
src/drt_xfem/enrichment_utils.cpp
src/drt_xfem/field_enriched.cpp
src/drt_xfem/interface.cpp
src/drt_xfem/interfacexfsi.cpp
src/drt_xfem/physics.cpp
src/drt_xfem/spacetime_boundary.cpp
src/drt_xfem/xdofmapcreation_combust.cpp
src/drt_xfem/xdofmapcreation.cpp
src/drt_xfem/xdofmapcreation_fsi.cpp
src/drt_xfem/xdofmapcreation_parallel_utils.cpp
)

set (input_full_source
src/input_full/inherit_design_dis.c
src/input_full/inherit_design_ele.c
src/input_full/inherit_insidedesign.c
src/input_full/inp_cal_dirich_fac.c
src/input_full/input_conditions.c
src/input_full/input_control_global.c
src/input_full/input_ctr_head.c
src/input_full/input_curves.c
src/input_full/input_design.c
src/input_full/input_design_top.c
src/input_full/input_funct.c
src/input_full/input_locsys.c
src/input_full/input_material.c
src/input_full/input_mesh.c
src/input_full/input_monitor.c
src/input_full/input_resultdescr.c
src/input_full/input_sol.c
src/input_full/input_submesh.c
src/input_full/input_surface_energy.c
src/input_full/input_topology.c
)

set (global_full_source
src/global_full/main_ccarat.c
src/global_full/global_control.c
src/global_full/global_inp_control.c
src/global_full/global_cal_control.c
src/global_full/global_init_control.c
src/global_full/global_ass_dof.c
src/global_full/global_ass_dof_ndis.c
src/global_full/global_dyn_control.c
src/global_full/global_timecurve.c
src/global_full/global_node_find.c
src/global_full/global_dofconnect.c
src/global_full/global_monitoring.c
src/global_full/global_locsys.c
src/global_full/global_result_test.c
src/global_full/global_calelm.c
src/global_full/global_calelm_fast.c
src/global_full/global_calrhs.c
src/global_full/global_restart_control.c
src/global_full/global_open_files.c
src/global_full/global_check_max.c
src/global_full/global_divide_fast.c
src/global_full/global_subdivide.c
src/global_full/global_subdivide_process.c
)

set (ccarat_math_source
src/math/arcline.c
src/math/geometry.c
src/math/math1.c
src/math/sort_find.c
)

set (ccarat_output_source
src/output/out_checkfilesize.c
src/output/out_gid_msh.c
src/output/out_gid_sol.c
src/output/out_gid_soldyn.c
src/output/out_gid_solfsi.c
src/output/out_gid_solssi.c
src/output/out_global.c
src/output/out_monitor.c
src/output/out_plt.c
)

set (ccarat_solver_source
src/solver/solver_add_data.c
src/solver/solver_add_data_fast.c
src/solver/solver_add_mlib.c
src/solver/solver_aztec.c
src/solver/solver_colsol.c
src/solver/solver_control.c
src/solver/solver_copy_matrix.c
src/solver/solver_hypre.c
src/solver/solver_lapack.c
src/solver/solver_mask_matrices.c
src/solver/solver_mask_mds.c
src/solver/solver_mask_service.c
src/solver/solver_mlib.c
src/solver/solver_mlpcg2.c
src/solver/solver_mlpcg3.c
src/solver/solver_mlpcg4.c
src/solver/solver_mlpcg_agg.c
src/solver/solver_mlpcg.c
src/solver/solver_mlpcg_csr.c
src/solver/solver_mlpcg_fish.c
src/solver/solver_mlpcg_matmat.c
src/solver/solver_mlpcg_overlap.c
src/solver/solver_mlpcg_smo.c
src/solver/solver_mlpcg_spooles.c
src/solver/solver_mlpcg_test.c
src/solver/solver_mumps.c
src/solver/solver_oll_add.c
src/solver/solver_oll_aztec.c
src/solver/solver_oll.c
src/solver/solver_oll_service2.c
src/solver/solver_oll_service.c
src/solver/solver_oll_spooles.c
src/solver/solver_service2.c
src/solver/solver_service3.c
src/solver/solver_service.c
src/solver/solver_sparse_bdcsr.c
src/solver/solver_sparse.c
src/solver/solver_sparse_ccf.c
src/solver/solver_sparse_dense.c
src/solver/solver_sparse_msr.c
src/solver/solver_sparse_parcsr.c
src/solver/solver_sparse_rcptr.c
src/solver/solver_sparse_skyline.c
src/solver/solver_sparse_spooles.c
src/solver/solver_sparse_ucchb.c
src/solver/solver_spooles.c
src/solver/solver_superlu.c
src/solver/solver_umfpack.c
)

set (pss_full_source
src/pss_full/pss_am.c
src/pss_full/pss_ds.c
src/pss_full/pss_fr.c
src/pss_full/pss_pss.c
src/pss_full/pss_perf.c
src/pss_full/pss_visual.c
src/pss_full/pss_table.c
src/pss_full/pss_table_iter.c
src/pss_full/pss_set.c
src/pss_full/pss_parser.c
)

set (baci_source
${adapter_source}
${ale_source}
${bele3_source}
${combust_element_source}
${combust_source}
${constraint_element_source}
${constraint_source}
${contact_source}
${contactnew_source}
${elch_source}
${fem_general_source}
${fluid_element_source}
${fluid_source}
${fsi_source}
${geometry_source}
${inpar_source}
${inv_analysis_source}
${io_source}
${lib_element_source}
${lib_source}
${loma_source}
${mat_source}
${matelast_source}
${mortar_source}
${nurbs_discret_source}
${potential_source}
${scatra_source}
${statmech_source}
${structure_source}
${stru_multi_source}
${surfstress_source}
${thermo_main_source}
${thermo_timint_source}
${tsi_source}
${xdiff3_source}
${xfem_source}
${input_full_source}
${global_full_source}
${ccarat_math_source}
${ccarat_output_source}
${ccarat_solver_source}
${pss_full_source}
)

#   all possible elements

option (D_ALE         "ale element (2D+3D)" OFF)
if (D_ALE)
  add_definitions("-DD_ALE")
  set (baci_source
    ${baci_source}
    ${ale2_source}
    ${ale3_source}
    )
endif(D_ALE)

option (D_FLUID2      "fluid2" OFF)
if (D_FLUID2)
  add_definitions("-DD_FLUID2")
  set( f2_source
    src/drt_f2/fluid2.cpp
    src/drt_f2/fluid2_evaluate.cpp
    src/drt_f2/fluid2_genalpha_resVMM.cpp
    src/drt_f2/fluid2_impl.cpp
    src/drt_f2/fluid2_input.cpp
    src/drt_f2/fluid2_line.cpp
    src/drt_f2/fluid2_line_evaluate.cpp
    src/drt_f2/fluid2_nurbs.cpp
    src/drt_f2/fluid2_stationary.cpp
    src/drt_f2/fluid2_th.cpp
    src/drt_f2/fluid2_weak_dbc.cpp
    )
  set (baci_source
    ${baci_source}
    ${f2_source}
    )
endif(D_FLUID2)

option (D_FLUID3      "fluid3" OFF)
if (D_FLUID3)
  add_definitions("-DD_FLUID3")
  set( f3_source
    src/drt_f3/fluid3.cpp
    src/drt_f3/fluid3_evaluate.cpp
    src/drt_f3/fluid3_genalpha_resVMM.cpp
    src/drt_f3/fluid3_impl.cpp
    src/drt_f3/fluid3_input.cpp
    src/drt_f3/fluid3_line.cpp
    src/drt_f3/fluid3_line_evaluate.cpp
    src/drt_f3/fluid3_lin_impl.cpp
    src/drt_f3/fluid3_nurbs.cpp
    src/drt_f3/fluid3_stabilization.cpp
    src/drt_f3/fluid3_stationary.cpp
    src/drt_f3/fluid3_surface.cpp
    src/drt_f3/fluid3_surface_evaluate.cpp
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
    )
  set (baci_source
    ${baci_source}
    ${f3_source}
    )
endif(D_FLUID3)

option (D_SHELL8      "shell8" OFF)
if (D_SHELL8)
  add_definitions("-DD_SHELL8")
  set( s8_source
    src/drt_s8/shell8.cpp
    src/drt_s8/shell8_evaluate.cpp
    src/drt_s8/shell8_input.cpp
    src/drt_s8/shell8_line.cpp
    src/drt_s8/shell8_line_evaluate.cpp
    )
  set (baci_source
    ${baci_source}
    ${s8_source}
    )
endif(D_SHELL8)

option (D_WALL1       "wall1" OFF)
if (D_WALL1)
  add_definitions("-DD_WALL1")
  set( w1_source
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
  set (baci_source
    ${baci_source}
    ${w1_source}
    )
endif(D_WALL1)

option (D_THERMO      "thermal elements" OFF)
if (D_THERMO)
  add_definitions("-DD_THERMO")
  set( thermo_element_source
    src/drt_thermo/thermo_ele_boundary_impl.cpp
    src/drt_thermo/thermo_ele_impl.cpp
    src/drt_thermo/thermo_ele_impl_utils.cpp
    src/drt_thermo/thermo_element_boundary_evaluate.cpp
    src/drt_thermo/thermo_element.cpp
    src/drt_thermo/thermo_element_evaluate.cpp
    src/drt_thermo/thermo_element_input.cpp
    )
  set (baci_source
    ${baci_source}
    ${thermo_element_source}
    )
endif(D_THERMO)

option (D_SOLID3      "spatial solid elements" OFF)
if (D_SOLID3)
  add_definitions("-DD_SOLID3")
  set( so3_source
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
    src/drt_so3/so_hex8p1j1.cpp
    src/drt_so3/so_hex8p1j1_evaluate.cpp
    src/drt_so3/so_hex8p1j1_input.cpp
    src/drt_so3/so_hex8_service.cpp
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
    src/drt_so3/so_sh8p8_auxiliar.cpp
    src/drt_so3/so_sh8p8.cpp
    src/drt_so3/so_sh8p8_eas.cpp
    src/drt_so3/so_sh8p8_evaluate.cpp
    src/drt_so3/so_sh8p8_input.cpp
    src/drt_so3/so_sh8_service.cpp
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
  set (baci_source
    ${baci_source}
    ${so3_source}
    )
endif(D_SOLID3)

option (D_BEAM2       "2D Timoshenko beam element" OFF)
if (D_BEAM2)
  add_definitions("-DD_BEAM2")
  set( beam2_source
    src/drt_beam2/beam2.cpp
    src/drt_beam2/beam2_evaluate.cpp
    src/drt_beam2/beam2_input.cpp
    )
  set (baci_source
    ${baci_source}
    ${beam2_source}
    )
endif(D_BEAM2)

option (D_BEAM2R      "2D Reissner beam element" OFF)
if (D_BEAM2R)
  add_definitions("-DD_BEAM2R")
  set( beam2r_source
    src/drt_beam2r/beam2r.cpp
    src/drt_beam2r/beam2r_evaluate.cpp
    src/drt_beam2r/beam2r_input.cpp
    )
  set (baci_source
    ${baci_source}
    ${beam2r_source}
    )
endif(D_BEAM2R)

option (D_BEAM3       "3D Reissner beam element" OFF)
if (D_BEAM3)
  add_definitions("-DD_BEAM3")
  set( beam3_source
    src/drt_beam3/beam3.cpp
    src/drt_beam3/beam3_evaluate.cpp
    src/drt_beam3/beam3_input.cpp
    )
  set (baci_source
    ${baci_source}
    ${beam3_source}
    )
endif(D_BEAM3)

option (D_TRUSS3      "3D truss element" OFF)
if (D_TRUSS3)
  add_definitions("-DD_TRUSS3")
  set( truss3_source
    src/drt_truss3/truss3.cpp
    src/drt_truss3/truss3_evaluate.cpp
    src/drt_truss3/truss3_input.cpp
    )
  set (baci_source
    ${baci_source}
    ${truss3_source}
    )
endif(D_TRUSS3)

option (D_TRUSS2      "2D truss element" OFF)
if (D_TRUSS2)
  add_definitions("-DD_TRUSS2")
  set( truss2_source
    src/drt_truss2/truss2.cpp
    src/drt_truss2/truss2_evaluate.cpp
    src/drt_truss2/truss2_input.cpp
    )
  set (baci_source
    ${baci_source}
    ${truss2_source}
    )
endif(D_TRUSS2)

option (D_TORSION3    "3D torsion spring element" OFF)
if (D_TORSION3)
  add_definitions("-DD_TORSION3")
  set( torsion3_source
    src/drt_torsion3/torsion3.cpp
    src/drt_torsion3/torsion3_evaluate.cpp
    src/drt_torsion3/torsion3_input.cpp
    )
  set (baci_source
    ${baci_source}
    ${torsion3_source}
    )
endif(D_TORSION3)

option (D_TORSION2    "2D torsion spring element" OFF)
if (D_TORSION2)
  add_definitions("-DD_TORSION2")
  set( torsion2_source
    src/drt_torsion2/torsion2.cpp
    src/drt_torsion2/torsion2_evaluate.cpp
    src/drt_torsion2/torsion2_input.cpp
    )
  set (baci_source
    ${baci_source}
    ${torsion2_source}
    )
endif(D_TORSION2)

option (D_ARTNET      "1D Arterial network elements" OFF)
if (D_ARTNET)
  add_definitions("-DD_ARTNET")
  set( art_net_source
    src/drt_art_net/artery.cpp
    src/drt_art_net/artery_evaluate.cpp
    src/drt_art_net/artery_input.cpp
    src/drt_art_net/artery_lin_exp.cpp
    src/drt_art_net/art_junction.cpp
    src/drt_art_net/art_net_dyn_drt.cpp
    src/drt_art_net/artnetexplicitintegration.cpp
    src/drt_art_net/art_terminal_bc.cpp
    src/drt_art_net/art_write_gnuplot.cpp
    )
  set (baci_source
    ${baci_source}
    ${art_net_source}
    )
endif(D_ARTNET)

option (D_RED_AIRWAYS "reduced dimensional airways for calculating flow in lungs" OFF)
if (D_RED_AIRWAYS)
  add_definitions("-DD_RED_AIRWAYS")
  set( red_airways_source
    src/drt_red_airways/airway_evaluate.cpp
    src/drt_red_airways/airway_impl.cpp
    src/drt_red_airways/airwayimplicitintegration.cpp
    src/drt_red_airways/red_airway.cpp
    src/drt_red_airways/red_airway_input.cpp
    src/drt_red_airways/red_airways_dyn_drt.cpp
    )
  set (baci_source
    ${baci_source}
    ${red_airways_source}
    )
endif(D_RED_AIRWAYS)

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

option (QHULL "needed for XFEM simulations" OFF)
if (QHULL)
  add_definitions("-DQHULL")
endif(QHULL)

option (D_EXODUS "Pre_exodus (exodusII, netcdf)" OFF)
if (D_EXODUS)
  add_definitions("-DD_EXODUS")
endif(D_EXODUS)

option (PRESTRESS "prestress a structure with Neumann loads only" OFF)
if (PRESTRESS)
  add_definitions("-DPRESTRESS")
endif(PRESTRESS)

option (POSTSTRESS "restart with a prestressed structure" OFF)
if (POSTSTRESS)
  add_definitions("-DPOSTSTRESS")
endif(POSTSTRESS)

option (INVERSEDESIGNCREATE "run an inverse design problem" OFF)
if (INVERSEDESIGNCREATE)
  add_definitions("-DINVERSEDESIGNCREATE")
endif(INVERSEDESIGNCREATE)

option (INVERSEDESIGNUSE "use a previously created inverse design" OFF)
if (INVERSEDESIGNUSE)
  add_definitions("-DINVERSEDESIGNUSE")
endif(INVERSEDESIGNUSE)

option (PARMETIS "use parmetis for parallel setup of node graph" OFF)
if (PARMETIS)
  add_definitions("-DPARMETIS")
endif(PARMETIS)

option (D_ALE_BFLOW "adds the effect of ALE defformation on the volumetric flow rate calculation on the boundaries" OFF)
if (D_ALE_BFLOW)
  add_definitions("-DD_ALE_BFLOW")
endif(D_ALE_BFLOW)


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

option (TRAP_FE "crash ccarat if a nan of inf occurs" OFF)
if (TRAP_FE)
  add_definitions("-DTRAP_FE")
endif(TRAP_FE)

option (THROWELEMENTERRORS "throw errors at element level so all element error can be found" OFF)
if (THROWELEMENTERRORS)
  add_definitions("-DTHROWELEMENTERRORS")
endif(THROWELEMENTERRORS)


#

add_executable (baci ${baci_source})

#find_path( BLITZ_CONFIG_DIR blitz/gnu/bzconfig.h
#  /usr/lib/blitz/include
#)

#find_path( NETCDF_INCLUDE_DIR netcdf.h
#   /usr/include/netcdf
#   /usr/include/netcdf-3
#)

#find_path( NETCDF_LIBRARY libnetcdf.so
#   /usr/lib
#)

#find_path( EXODUS_INCLUDE_DIR exodusII.h
#  /lnm/lib/fc8/openmpi/include
#)

#find_path( EXODUS_LIBRARY libexoIIv2c.a
#  /lnm/lib/fc8/openmpi/lib
#)

#find_path( NEMESIS_INCLUDE_DIR ne_nemesisI.h
#  /lnm/lib/fc8/openmpi/include
#)

#find_path( NEMESIS_LIBRARY libnemIc.a
#  /lnm/lib/fc8/openmpi/lib
#)

#find_path( METIS_LIBRARY libmetis.a
#  /lnm/lib/fc8/openmpi/lib
#)

#find_path( PARMETIS_LIBRARY libparmetis.a
#  /lnm/lib/fc8/openmpi/lib
#)

#message(${LIBRARIES})

target_link_libraries (baci
  ${LIBRARIES}
#global_full
##
## multi field algorithms
#drt_fsi
#drt_loma
#drt_elch
#drt_thermo_main
##
## adapter
#drt_adapter
##
## single field algorithms
#drt_scatra
#drt_ale
#drt_art_net
#drt_combust
#drt_fluid
#drt_io
#drt_nurbs_discret
#drt_structure
#drt_inv_analysis
#drt_statmech
#drt_surfstress
#drt_thermo_timint
#drt_tsi
##
## support
#drt_potential
#drt_contact
#drt_contactnew
#drt_mortar
##
## element creation part of framework
#drt_lib_element
##
## elements
#drt_torsion2
#drt_torsion3
#drt_truss2
#drt_truss3
#drt_f2
#drt_f3
#drt_s8
#drt_so3
#drt_beam2
#drt_beam2r
#drt_beam3
#drt_bele3
#drt_ale2
#drt_ale3
#drt_thermo_element
#drt_constraint_element
#drt_combust_element
#drt_w1
#drt_xdiff3
#drt_xfem
##
## support code needed by elements
#drt_fluid_element
##drt_stru_multi
##
## framework
#drt_constraint
#drt_geometry
#drt_fem_general
#drt_lib
#drt_mat
#drt_matelast
#drt_inpar
#pss_full
#ccarat_output
#ccarat_math
#ccarat_solver
#input_full
#
)
