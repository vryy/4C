/*----------------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problem with non-matching grids using a monolithic scheme
in saddle-point formulation with Lagrange multipliers discretized on the fluid interface

\level 1
*/

/*----------------------------------------------------------------------------*/


#include "4C_fsi_mortarmonolithic_fluidsplit_sp.hpp"

#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_constraint_manager.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_coupling_adapter_mortar.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_debugwriter.hpp"
#include "4C_fsi_statustest.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linalg_utils_sparse_algebra_print.hpp"
#include "4C_mortar_utils.hpp"
#include "4C_structure_aux.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
FSI::MortarMonolithicFluidSplitSaddlePoint::MortarMonolithicFluidSplitSaddlePoint(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : BlockMonolithic(comm, timeparams), comm_(comm)
{
  // ---------------------------------------------------------------------------
  // FSI specific check of Dirichlet boundary conditions
  // ---------------------------------------------------------------------------
  // Create intersection of slave DOFs that hold a Dirichlet boundary condition
  // and are located at the FSI interface
  std::vector<Teuchos::RCP<const Epetra_Map>> intersectionmaps;
  intersectionmaps.push_back(fluid_field()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(fluid_field()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      Core::LinAlg::MultiMapExtractor::IntersectMaps(intersectionmaps);

  // Check whether the intersection is empty
  if (intersectionmap->NumGlobalElements() != 0)
  {
    //      std::cout << "Slave interface nodes with Dirichlet boundary condition "
    //                "(input file numbering):" << std::endl;
    //      for (int i=0; i < (int)fluid_field()->discretization()->NumMyRowNodes(); i++)
    //      {
    //        // get all nodes and add them
    //        int gid = fluid_field()->discretization()->NodeRowMap()->GID(i);

    //        // do only nodes that I have in my discretization
    //        if (!fluid_field()->discretization()->NodeColMap()->MyGID(gid)) continue;
    //        Core::Nodes::Node* node = fluid_field()->discretization()->gNode(gid);
    //        if (!node) FOUR_C_THROW("Cannot find node with gid %",gid);

    //        std::vector<int> nodedofs = fluid_field()->discretization()->Dof(node);

    //        for (int j=0; j < (int)nodedofs.size(); j++)
    //        {
    //          for (int k=0; k < (int)intersectionmap->NumGlobalElements(); k++)
    //          {
    //            if (nodedofs[j] == intersectionmap->GID(k))
    //            {
    //              std::cout << gid+1 << std::endl;
    //              k = (int)intersectionmap->GID(k);
    //              j = (int)nodedofs.size();
    //            }
    //          }
    //        }
    //      }

    // It is not allowed, that slave DOFs at the interface hold a Dirichlet
    // boundary condition. Thus --> Error message

    // We do not have to care whether ALE interface DOFs carry DBCs in the
    // input file since they do not occur in the monolithic system and, hence,
    // do not cause a conflict.

    std::stringstream errormsg;
    errormsg << "  "
                "+---------------------------------------------------------------------------------"
                "------------+"
             << std::endl
             << "  |                DIRICHLET BOUNDARY CONDITIONS ON SLAVE SIDE OF FSI INTERFACE   "
                "              |"
             << std::endl
             << "  "
                "+---------------------------------------------------------------------------------"
                "------------+"
             << std::endl
             << "  | NOTE: The slave side of the interface is not allowed to carry Dirichlet "
                "boundary conditions.|"
             << std::endl
             << "  |                                                                               "
                "              |"
             << std::endl
             << "  | This is a fluid split scheme. Hence, master and slave field are chosen as "
                "follows:          |"
             << std::endl
             << "  |     MASTER  = STRUCTURE                                                       "
                "              |"
             << std::endl
             << "  |     SLAVE   = FLUID                                                           "
                "              |"
             << std::endl
             << "  |                                                                               "
                "              |"
             << std::endl
             << "  | Dirichlet boundary conditions were detected on slave interface degrees of "
                "freedom. Please   |"
             << std::endl
             << "  | remove Dirichlet boundary conditions from the slave side of the FSI "
                "interface.              |"
             << std::endl
             << "  | Only the master side of the FSI interface is allowed to carry Dirichlet "
                "boundary conditions.|"
             << std::endl
             << "  "
                "+---------------------------------------------------------------------------------"
                "------------+"
             << std::endl;

    FOUR_C_THROW(errormsg.str());
  }
  // ---------------------------------------------------------------------------

  notsetup_ = true;

  coupling_solid_fluid_mortar_ = Teuchos::rcp(new Core::Adapter::CouplingMortar(
      Global::Problem::Instance()->NDim(), Global::Problem::Instance()->mortar_coupling_params(),
      Global::Problem::Instance()->contact_dynamic_params(),
      Global::Problem::Instance()->spatial_approximation_type()));

  ale_inner_interf_transform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);
  fluid_mesh_inner_inner_transform_ = Teuchos::rcp(new Core::LinAlg::MatrixColTransform);

  create_lagrange_multiplier_dof_row_map();
  set_lag_mult();

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (coupling_solid_fluid_mortar_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'coupling_solid_fluid_mortar_' failed.");
  }
  if (ale_inner_interf_transform_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'ale_inner_interf_transform_' failed.");
  }
  if (fluid_mesh_inner_inner_transform_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'fluid_mesh_inner_inner_transform_' failed.");
  }
  if (lag_mult_ == Teuchos::null)
  {
    FOUR_C_THROW("Allocation of 'lag_mult_' failed.");
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::set_lag_mult()
{
  lag_mult_ = Teuchos::rcp(new Epetra_Vector(*lag_mult_dof_map_, true));
  lag_mult_old_ = Teuchos::rcp(new Epetra_Vector(*lag_mult_dof_map_, true));
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::SetupSystem()
{
  if (notsetup_)
  {
    const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();

    set_default_parameters(fsidyn, nox_parameter_list());

    // we use non-matching meshes at the interface
    // mortar with: structure = master, fluid = slave

    const int ndim = Global::Problem::Instance()->NDim();

    // get coupling objects
    Core::Adapter::Coupling& interface_coup_fluid_ale = interface_fluid_ale_coupling();

    /* structure to fluid
     * coupling condition at the fsi interface:
     * displacements (=number spacial dimensions) are coupled
     * e.g.: 3D: coupleddof = [1, 1, 1]
     */
    std::vector<int> coupleddof(ndim, 1);

    coupling_solid_fluid_mortar_->Setup(structure_field()->discretization(),
        fluid_field()->discretization(), ale_field()->write_access_discretization(), coupleddof,
        "FSICoupling", comm_, Global::Problem::Instance()->FunctionManager(), true);

    // fluid to ale at the interface
    interface_coup_fluid_ale.setup_condition_coupling(*fluid_field()->discretization(),
        fluid_field()->Interface()->FSICondMap(), *ale_field()->discretization(),
        ale_field()->Interface()->FSICondMap(), "FSICoupling", ndim);

    Core::Adapter::Coupling& coup_fluid_ale = fluid_ale_coupling();

    // the fluid-ale coupling always matches
    const Epetra_Map* fluidnodemap = fluid_field()->discretization()->NodeRowMap();
    const Epetra_Map* alenodemap = ale_field()->discretization()->NodeRowMap();

    coup_fluid_ale.setup_coupling(*fluid_field()->discretization(), *ale_field()->discretization(),
        *fluidnodemap, *alenodemap, ndim);

    fluid_field()->SetMeshMap(coup_fluid_ale.MasterDofMap());

    create_combined_dof_row_map();

    /*------------------------------------------------------------------------*/
    // Switch fluid to interface split block matrix
    fluid_field()->use_block_matrix(true);

    // build ale system matrix in splitted system
    ale_field()->create_system_matrix(ale_field()->Interface());

    aleresidual_ = Teuchos::rcp(new Epetra_Vector(*ale_field()->Interface()->OtherMap()));

    // -------------------------------------------------------------------------
    // Build the global Dirichlet map extractor
    setup_dbc_map_extractor();
    // -------------------------------------------------------------------------#

    // enable debugging
    if (Core::UTILS::IntegralValue<int>(fsidyn, "DEBUGOUTPUT") & 2)
    {
      pcdbg_ = Teuchos::rcp(new UTILS::MonolithicDebugWriter(*this));
    }

    create_system_matrix();

    notsetup_ = false;
  }

  const int restart = Global::Problem::Instance()->Restart();
  if (restart)
  {
    const bool restartfrompartfsi =
        Core::UTILS::IntegralValue<bool>(timeparams_, "RESTART_FROM_PART_FSI");
    if (restartfrompartfsi)  // restart from part. fsi
    {
      if (comm_.MyPID() == 0)
        std::cout << "Warning: RESTART_FROM_PART_FSI for mortar fsi is not jet implemented. For "
                     "now lambda_ is simply assumed to be zero!"
                  << std::endl;

      // mortar business still has to be done here..
    }
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::create_lagrange_multiplier_dof_row_map()
{
  const int num_glob_elem_fluid_interface =
      fluid_field()->Interface()->FSICondMap()->NumGlobalElements();
  const int num_loc_elem_fluid_interface =
      fluid_field()->Interface()->FSICondMap()->NumMyElements();
  const int max_gid_ale = ale_field()->dof_row_map()->MaxAllGID();
  lag_mult_dof_map_ = Teuchos::rcp(new Epetra_Map(
      num_glob_elem_fluid_interface, num_loc_elem_fluid_interface, max_gid_ale + 1, comm_));
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::create_combined_dof_row_map()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(structure_field()->dof_row_map());
  vecSpaces.push_back(fluid_field()->dof_row_map());
  vecSpaces.push_back(ale_field()->Interface()->OtherMap());
  vecSpaces.push_back(lag_mult_dof_map_);

  if (vecSpaces[1]->NumGlobalElements() == 0)
    FOUR_C_THROW("No inner fluid equations. Splitting not possible.");

  set_dof_row_maps(vecSpaces);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::create_system_matrix()
{
  FSI::BlockMonolithic::create_system_matrix(systemmatrix_, false);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Combo>
FSI::MortarMonolithicFluidSplitSaddlePoint::create_status_test(
    Teuchos::ParameterList& nlParams, Teuchos::RCP<::NOX::Epetra::Group> grp)
{
  // ---------------------------------------------------------------------------
  // Setup the test framework
  // ---------------------------------------------------------------------------
  // Create the top-level test combo
  Teuchos::RCP<::NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::OR));

  // Create test combo for convergence of residuals and iterative increments
  Teuchos::RCP<::NOX::StatusTest::Combo> converged =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::AND));

  // Create some other plausibility tests
  Teuchos::RCP<::NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new ::NOX::StatusTest::MaxIters(nlParams.get<int>("Max Iterations")));
  Teuchos::RCP<::NOX::StatusTest::FiniteValue> fv =
      Teuchos::rcp(new ::NOX::StatusTest::FiniteValue);

  // Add single tests to the top-level test combo
  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // Start filling the 'converged' combo here
  // require one solve
  converged->addStatusTest(Teuchos::rcp(new NOX::FSI::MinIters(1)));

  // ---------------------------------------------------------------------------
  // setup tests for structural displacement field
  // ---------------------------------------------------------------------------
  // create ::NOX::StatusTest::Combo for structural displacement field
  Teuchos::RCP<::NOX::StatusTest::Combo> structcombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp_L2 = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "DISPL residual", extractor(), 0, nlParams.get<double>("Tol dis res L2"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "DISPL residual", extractor(), 0, nlParams.get<double>("Tol dis res Inf"),
      ::NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("DISPL update", extractor(), 0,
          nlParams.get<double>("Tol dis inc L2"), ::NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("DISPL update", extractor(), 0,
          nlParams.get<double>("Tol dis inc Inf"), ::NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  add_status_test(structureDisp_L2);

  // add norm-tests to structural displacement ::NOX::StatusTest::Combo
  structcombo->addStatusTest(structureDisp_L2);
  structcombo->addStatusTest(structureDisp_inf);
  structcombo->addStatusTest(structureDispUpdate_L2);
  structcombo->addStatusTest(structureDispUpdate_inf);

  // add structural displacement test combo to top-level test combo
  converged->addStatusTest(structcombo);
  // ---------- end of structural displacement field tests

  // ---------------------------------------------------------------------------
  // setup tests for fluid velocity field
  // ---------------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map>> fluidvel;
  fluidvel.push_back(fluid_field()->InnerVelocityRowMap());
  fluidvel.push_back(Teuchos::null);
  Core::LinAlg::MultiMapExtractor fluidvelextract(*dof_row_map(), fluidvel);

  // create ::NOX::StatusTest::Combo for fluid velocity field
  Teuchos::RCP<::NOX::StatusTest::Combo> fluidvelcombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel_L2 = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "VELOC residual", fluidvelextract, 0, nlParams.get<double>("Tol vel res L2"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "VELOC residual", fluidvelextract, 0, nlParams.get<double>("Tol vel res Inf"),
      ::NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("VELOC update", fluidvelextract, 0,
          nlParams.get<double>("Tol vel inc L2"), ::NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("VELOC update", fluidvelextract, 0,
          nlParams.get<double>("Tol vel inc Inf"), ::NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  add_status_test(innerFluidVel_L2);

  // add norm-tests to fluid velocity ::NOX::StatusTest::Combo
  fluidvelcombo->addStatusTest(innerFluidVel_L2);
  fluidvelcombo->addStatusTest(innerFluidVel_inf);
  fluidvelcombo->addStatusTest(innerFluidVelUpdate_L2);
  fluidvelcombo->addStatusTest(innerFluidVelUpdate_inf);

  // add fluid velocity test combo to top-level test combo
  converged->addStatusTest(fluidvelcombo);
  // ---------- end of fluid velocity field tests

  // ---------------------------------------------------------------------------
  // setup tests for fluid pressure field
  // ---------------------------------------------------------------------------
  // build mapextractor
  std::vector<Teuchos::RCP<const Epetra_Map>> fluidpress;
  fluidpress.push_back(fluid_field()->PressureRowMap());
  fluidpress.push_back(Teuchos::null);
  Core::LinAlg::MultiMapExtractor fluidpressextract(*dof_row_map(), fluidpress);

  // create ::NOX::StatusTest::Combo for fluid pressure field
  Teuchos::RCP<::NOX::StatusTest::Combo> fluidpresscombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress_L2 = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "PRESS residual", fluidpressextract, 0, nlParams.get<double>("Tol pre res L2"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "PRESS residual", fluidpressextract, 0, nlParams.get<double>("Tol pre res Inf"),
      ::NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("PRESS update", fluidpressextract, 0,
          nlParams.get<double>("Tol pre inc L2"), ::NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("PRESS update", fluidpressextract, 0,
          nlParams.get<double>("Tol pre inc Inf"), ::NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  add_status_test(fluidPress_L2);

  // add norm-tests to fluid pressure ::NOX::StatusTest::Combo
  fluidpresscombo->addStatusTest(fluidPress_L2);
  fluidpresscombo->addStatusTest(fluidPress_inf);
  fluidpresscombo->addStatusTest(fluidPressUpdate_L2);
  fluidpresscombo->addStatusTest(fluidPressUpdate_inf);

  // add fluid pressure test combo to top-level test combo
  converged->addStatusTest(fluidpresscombo);
  // ---------- end of fluid pressure field tests

  // ---------------------------------------------------------------------------
  // setup tests for Lagrange multiplier field
  // ---------------------------------------------------------------------------
  // create ::NOX::StatusTest::Combo for Lagrange multiplier field
  Teuchos::RCP<::NOX::StatusTest::Combo> lag_mult_combo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::AND));

  // create Norm-objects for each norm that has to be tested
  Teuchos::RCP<NOX::FSI::PartialNormF> lag_mult_L2 = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "LAGMULT residual", extractor(), 3, nlParams.get<double>("Tol fsi res L2"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> lag_mult_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "LAGMULT residual", extractor(), 3, nlParams.get<double>("Tol fsi res Inf"),
      ::NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> lag_mult_update_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("LAGMULT update", extractor(), 3,
          nlParams.get<double>("Tol fsi inc L2"), ::NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> lag_mult_update_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("LAGMULT update", extractor(), 3,
          nlParams.get<double>("Tol fsi inc Inf"), ::NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  add_status_test(lag_mult_L2);

  // add norm-tests to structural displacement ::NOX::StatusTest::Combo
  lag_mult_combo->addStatusTest(lag_mult_L2);
  lag_mult_combo->addStatusTest(lag_mult_inf);
  lag_mult_combo->addStatusTest(lag_mult_update_L2);
  lag_mult_combo->addStatusTest(lag_mult_update_inf);

  // add structural displacement test combo to top-level test combo
  converged->addStatusTest(lag_mult_combo);
  // ---------- end of Lagrange multiplier field tests

  // Finally, return the test combo
  return combo;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::setup_dbc_map_extractor()
{
  /* Dirichlet maps for structure and fluid do not intersect with interface map.
   * ALE Dirichlet map might intersect with interface map, but ALE interface
   * DOFs are not part of the final system of equations. Hence, we just need the
   * intersection of inner ALE DOFs with Dirichlet ALE DOFs.
   */
  std::vector<Teuchos::RCP<const Epetra_Map>> aleintersectionmaps;
  aleintersectionmaps.push_back(ale_field()->GetDBCMapExtractor()->CondMap());
  aleintersectionmaps.push_back(ale_field()->Interface()->OtherMap());
  Teuchos::RCP<const Epetra_Map> aleintersectionmap =
      Core::LinAlg::MultiMapExtractor::IntersectMaps(aleintersectionmaps);

  // Merge Dirichlet maps of structure, fluid and ALE to global FSI Dirichlet map
  std::vector<Teuchos::RCP<const Epetra_Map>> dbcmaps;
  dbcmaps.push_back(structure_field()->GetDBCMapExtractor()->CondMap());
  dbcmaps.push_back(fluid_field()->GetDBCMapExtractor()->CondMap());
  dbcmaps.push_back(aleintersectionmap);

  Teuchos::RCP<const Epetra_Map> dbcmap = Core::LinAlg::MultiMapExtractor::MergeMaps(dbcmaps);

  // Finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new Core::LinAlg::MapExtractor(*dof_row_map(), dbcmap, true));
  if (dbcmaps_ == Teuchos::null) FOUR_C_THROW("Creation of FSI Dirichlet map extractor failed.");
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>
FSI::MortarMonolithicFluidSplitSaddlePoint::SystemMatrix() const
{
  return systemmatrix_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::initial_guess(
    Teuchos::RCP<Epetra_Vector> initial_guess)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MortarMonolithicFluidSplitSaddlePoint::initial_guess");

  Teuchos::RCP<const Epetra_Vector> lag_mult_initial_guess =
      Teuchos::rcp(new Epetra_Vector(*lag_mult_dof_map_, true));

  combine_field_vectors(*initial_guess, structure_field()->initial_guess(),
      fluid_field()->initial_guess(), ale_field()->initial_guess(), lag_mult_initial_guess, true);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::combine_field_vectors(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> solid_vector, Teuchos::RCP<const Epetra_Vector> fluid_vector,
    Teuchos::RCP<const Epetra_Vector> ale_vector, Teuchos::RCP<const Epetra_Vector> lag_mult_vector,
    bool fullvectors)
{
  if (fullvectors)
  {
    // extract inner DOFs from slave vectors
    Teuchos::RCP<const Epetra_Vector> ale_other_vector =
        ale_field()->Interface()->ExtractOtherVector(ale_vector);

    extractor().AddVector(*solid_vector, 0, f);
    extractor().AddVector(*fluid_vector, 1, f);
    extractor().AddVector(*ale_other_vector, 2, f);
    extractor().AddVector(*lag_mult_vector, 3, f);
  }
  else
  {
    extractor().AddVector(*solid_vector, 0, f);
    extractor().AddVector(*fluid_vector, 1, f);
    extractor().AddVector(*ale_vector, 2, f);
    extractor().AddVector(*lag_mult_vector, 3, f);
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::setup_rhs_residual(Epetra_Vector& f)
{
  // get single field residuals
  Teuchos::RCP<const Epetra_Vector> solid_single_field_rhs_vector =
      Teuchos::rcp(new Epetra_Vector(*structure_field()->RHS()));
  Teuchos::RCP<const Epetra_Vector> fluid_single_field_rhs_vector =
      Teuchos::rcp(new Epetra_Vector(*fluid_field()->RHS()));
  Teuchos::RCP<const Epetra_Vector> ale_single_field_rhs_vector =
      Teuchos::rcp(new Epetra_Vector(*ale_field()->RHS()));
  Teuchos::RCP<Epetra_Vector> lag_mult_rhs_vector =
      Teuchos::rcp(new Epetra_Vector(*lag_mult_dof_map_, true));

  // put the single field residuals together
  combine_field_vectors(f, solid_single_field_rhs_vector, fluid_single_field_rhs_vector,
      ale_single_field_rhs_vector, lag_mult_rhs_vector, true);

  // add additional ale residual to avoid incremental ale errors
  extractor().AddVector(*aleresidual_, 2, f);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::setup_rhs_lambda(Epetra_Vector& f)
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double solid_time_int_param = structure_field()->TimIntParam();
  const double fluid_time_int_param = fluid_field()->TimIntParam();
  const double fluid_res_scale = fluid_field()->residual_scaling();

  // get the mortar structure to fluid coupling matrix M
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortar_m =
      coupling_solid_fluid_mortar_->GetMortarMatrixM();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mortar_m_transf =
      Mortar::MatrixRowTransformGIDs(mortar_m, lag_mult_dof_map_);

  // get the mortar fluid to structure coupling matrix D
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortar_d =
      coupling_solid_fluid_mortar_->GetMortarMatrixD();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mortar_d_transf =
      Mortar::MatrixRowTransformGIDs(mortar_d, lag_mult_dof_map_);

  Teuchos::RCP<Epetra_Vector> lag_mult_step_increment =
      Teuchos::rcp(new Epetra_Vector(*lag_mult_dof_map_, true));
  lag_mult_step_increment->Update(1.0, *lag_mult_, -1.0, *lag_mult_old_, 0.0);

  // helper variables
  Teuchos::RCP<Epetra_Vector> lag_mult_old_rhs_struc_interf =
      Teuchos::rcp(new Epetra_Vector(mortar_m_transf->DomainMap(), true));
  Teuchos::RCP<Epetra_Vector> lag_mult_old_rhs_fluid_interf =
      Teuchos::rcp(new Epetra_Vector(mortar_d_transf->DomainMap(), true));

  mortar_m_transf->Multiply(true, *lag_mult_old_, *lag_mult_old_rhs_struc_interf);
  mortar_d_transf->Multiply(true, *lag_mult_old_, *lag_mult_old_rhs_fluid_interf);

  Teuchos::RCP<Epetra_Vector> lag_mult_old_rhs_struc_interf_full =
      structure_field()->Interface()->InsertFSICondVector(lag_mult_old_rhs_struc_interf);
  Teuchos::RCP<Epetra_Vector> lag_mult_old_rhs_fluid_interf_full =
      fluid_field()->Interface()->InsertFSICondVector(lag_mult_old_rhs_fluid_interf);

  lag_mult_old_rhs_fluid_interf_full->Scale(-1.0 / fluid_res_scale);

  // add lagrange multiplier
  extractor().AddVector(*lag_mult_old_rhs_struc_interf_full, 0, f);
  extractor().AddVector(*lag_mult_old_rhs_fluid_interf_full, 1, f);

  // helper variables
  Teuchos::RCP<Epetra_Vector> lag_mult_step_increment_rhs_struc_interf =
      Teuchos::rcp(new Epetra_Vector(mortar_m_transf->DomainMap(), true));
  Teuchos::RCP<Epetra_Vector> lag_mult_step_increment_rhs_fluid_interf =
      Teuchos::rcp(new Epetra_Vector(mortar_d_transf->DomainMap(), true));

  mortar_m_transf->Multiply(
      true, *lag_mult_step_increment, *lag_mult_step_increment_rhs_struc_interf);
  mortar_d_transf->Multiply(
      true, *lag_mult_step_increment, *lag_mult_step_increment_rhs_fluid_interf);

  Teuchos::RCP<Epetra_Vector> lag_mult_step_increment_rhs_struc_interf_full =
      structure_field()->Interface()->InsertFSICondVector(lag_mult_step_increment_rhs_struc_interf);
  Teuchos::RCP<Epetra_Vector> lag_mult_step_increment_rhs_fluid_interf_full =
      fluid_field()->Interface()->InsertFSICondVector(lag_mult_step_increment_rhs_fluid_interf);

  lag_mult_step_increment_rhs_struc_interf_full->Scale(1.0 * (1. - solid_time_int_param));
  lag_mult_step_increment_rhs_fluid_interf_full->Scale(
      -1.0 * (1. - fluid_time_int_param) / fluid_res_scale);

  // add lagrange multiplier
  extractor().AddVector(*lag_mult_step_increment_rhs_struc_interf_full, 0, f);
  extractor().AddVector(*lag_mult_step_increment_rhs_fluid_interf_full, 1, f);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::setup_rhs_firstiter(Epetra_Vector& f)
{
  // old interface velocity of fluid field
  const Teuchos::RCP<const Epetra_Vector> fluid_veln = fluid_field()->extract_interface_veln();

  // get the mortar structure to fluid coupling matrix M
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortar_m =
      coupling_solid_fluid_mortar_->GetMortarMatrixM();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mortar_m_transf =
      Mortar::MatrixRowTransformGIDs(mortar_m, lag_mult_dof_map_);

  // get the mortar fluid to structure coupling matrix D
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortar_d =
      coupling_solid_fluid_mortar_->GetMortarMatrixD();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mortar_d_transf =
      Mortar::MatrixRowTransformGIDs(mortar_d, lag_mult_dof_map_);

  // get fluid shape derivatives matrix
  const Teuchos::RCP<const Core::LinAlg::BlockSparseMatrixBase> fluid_shape_deriv =
      fluid_field()->ShapeDerivatives();

  // get ale matrix
  const Teuchos::RCP<const Core::LinAlg::BlockSparseMatrixBase> aleblock =
      ale_field()->BlockSystemMatrix();

  // extract ale submatrix
  const Core::LinAlg::SparseMatrix& ale_inner_interf = aleblock->Matrix(0, 1);

  // right hand side of single set of DOFs
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null;

  /* Different contributions/terms to the rhs are separated by the following
   * comment line */
  // ---------- inner fluid DOFs
  /* The following term is added to the inner fluid DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  dt * F^{G}_{I \Gamma} * u^{n}_{\Gamma}
   *
   */
  // ----------addressing term 1
  if (fluid_shape_deriv != Teuchos::null)
  {
    const Core::LinAlg::SparseMatrix& fluid_mesh_inner_interf = fluid_shape_deriv->Matrix(0, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(fluid_mesh_inner_interf.RangeMap(), true));

    fluid_mesh_inner_interf.Apply(*fluid_veln, *rhs);

    rhs->Scale(Dt());

    rhs = fluid_field()->Interface()->InsertOtherVector(rhs);

    extractor().AddVector(*rhs, 1, f);
  }
  // ----------end of term 1
  // ----------end of inner fluid DOFs

  // ---------- interface fluid DOFs
  /* The following term is added to the interface fluid DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  dt * F^{G}_{\Gamma \Gamma} * u^{n}_{\Gamma}
   *
   */
  // ----------addressing term 1
  if (fluid_shape_deriv != Teuchos::null)
  {
    const Core::LinAlg::SparseMatrix& fluid_mesh_interf_interf = fluid_shape_deriv->Matrix(1, 1);
    rhs = Teuchos::rcp(new Epetra_Vector(fluid_mesh_interf_interf.RangeMap(), true));

    fluid_mesh_interf_interf.Apply(*fluid_veln, *rhs);
    rhs->Scale(Dt());
    rhs = fluid_field()->Interface()->InsertFSICondVector(rhs);

    extractor().AddVector(*rhs, 1, f);
  }
  // ----------end of term 1
  // ----------end of interface fluid DOFs

  // ---------- inner ALE DOFs
  /* The following term is added to the inner ALE DOFs of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  dt * A_{I \Gamma} * u^{n}_{\Gamma}
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(ale_inner_interf.RangeMap(), true));
  ale_inner_interf.Apply(*fluid_veln, *rhs);
  rhs->Scale(-1. * Dt());

  extractor().AddVector(*rhs, 2, f);
  // ----------end of term 1
  // ----------end of inner ALE DOFs

  // ---------- lagrange multiplier
  /* The following term is added to the lagrange multiplier of right hand side:
   *
   * rhs_firstnewtonstep =
   *
   * (1)  + dt * D * u^{n}_{\Gamma}
   *
   * (2)  - M * \Delta d_{\Gamma,p}
   *
   */
  // ----------addressing term 1
  rhs = Teuchos::rcp(new Epetra_Vector(*lag_mult_dof_map_, true));

  mortar_d_transf->Apply(*fluid_veln, *rhs);
  rhs->Scale(Dt());

  extractor().AddVector(*rhs, 3, f);
  // ----------end of term 1

  // ----------addressing term 2
  rhs = Teuchos::rcp(new Epetra_Vector(*lag_mult_dof_map_, true));

  mortar_m_transf->Apply(*ddgpred_, *rhs);
  rhs->Scale(-1.);

  extractor().AddVector(*rhs, 3, f);
  // ----------end of term 2
  // ----------end of lagrange multiplier
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::setup_system_matrix(
    Core::LinAlg::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MortarMonolithicFluidSplitSaddlePoint::setup_system_matrix");

  // get the mortar structure to fluid coupling matrix M
  const Teuchos::RCP<const Core::LinAlg::SparseMatrix> mortar_m =
      coupling_solid_fluid_mortar_->GetMortarMatrixM();

  // get the mortar fluid to structure coupling matrix D
  const Teuchos::RCP<const Core::LinAlg::SparseMatrix> mortar_d =
      coupling_solid_fluid_mortar_->GetMortarMatrixD();

  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double solid_time_int_param = structure_field()->TimIntParam();
  const double fluid_time_int_param = fluid_field()->TimIntParam();
  const double fluid_res_scale = fluid_field()->residual_scaling();

  // time scaling factor for fluid
  const double fluid_timescale = fluid_field()->TimeScaling();

  // get fluid shape derivatives matrix
  const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> fluid_shape_deriv =
      fluid_field()->ShapeDerivatives();

  // get single field block matrices
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> solidblock = structure_field()->system_matrix();
  const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> fluidblock =
      fluid_field()->BlockSystemMatrix();
  const Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> aleblock =
      ale_field()->BlockSystemMatrix();

  // extract submatrices
  const Core::LinAlg::SparseMatrix& fluid_inner_inner = fluidblock->Matrix(0, 0);
  const Core::LinAlg::SparseMatrix& fluid_interf_inner = fluidblock->Matrix(1, 0);
  const Core::LinAlg::SparseMatrix& fluid_inner_interf = fluidblock->Matrix(0, 1);
  const Core::LinAlg::SparseMatrix& fluid_interf_interf = fluidblock->Matrix(1, 1);
  const Core::LinAlg::SparseMatrix& ale_inner_inner = aleblock->Matrix(0, 0);
  const Core::LinAlg::SparseMatrix& ale_inner_interf = aleblock->Matrix(0, 1);

  // ---------------------------------------------------------------------------
  // BEGIN building the global 6x6 system matrix
  // ---------------------------------------------------------------------------
  // Contributions to blocks in system matrix are listed separately.
  // Block numbering in comments ranges from (1,1) to (6,6).

  // ---------Addressing contribution to blocks (1,1),(1,2),(2,1),(2,2)
  mat.Assign(0, 0, Core::LinAlg::View, *solidblock);

  // ---------Addressing contribution to blocks (3,3),(3,4),(4,3),(4,4)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> aux_fluidblock =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(fluidblock->FullRowMap(), 108, false));
  aux_fluidblock->Add(fluid_inner_inner, false, 1.0, 0.0);
  aux_fluidblock->Add(fluid_interf_inner, false, 1.0, 1.0);
  aux_fluidblock->Add(fluid_inner_interf, false, 1.0, 1.0);
  aux_fluidblock->Add(fluid_interf_interf, false, 1.0, 1.0);
  aux_fluidblock->Complete(fluidblock->FullDomainMap(), fluidblock->FullRangeMap(), true);

  // ---------Addressing contribution to block (5,4)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> aux_ale_inner_interf =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(ale_inner_inner.RowMap(), 81, false));
  (*ale_inner_interf_transform_)(aleblock->FullRowMap(), aleblock->FullColMap(), ale_inner_interf,
      1., Core::Adapter::CouplingSlaveConverter(interface_fluid_ale_coupling()),
      *aux_ale_inner_interf);

  aux_ale_inner_interf->Scale(1. / fluid_timescale);
  aux_ale_inner_interf->Complete(fluidblock->DomainMap(), aux_ale_inner_interf->RangeMap(), true);

  mat.Assign(2, 1, Core::LinAlg::View, *aux_ale_inner_interf);

  // ---------Addressing contribution to block (5,5)
  mat.Assign(2, 2, Core::LinAlg::View, ale_inner_inner);

  // ---------Addressing contribution to block (6,2)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> aux_mortar_m =
      Mortar::MatrixRowTransformGIDs(mortar_m, lag_mult_dof_map_);
  aux_mortar_m->Complete(solidblock->DomainMap(), *lag_mult_dof_map_, true);

  mat.Assign(3, 0, Core::LinAlg::View, *aux_mortar_m);

  // ---------Addressing contribution to block (2,6)
  aux_mortar_m = Mortar::MatrixRowTransformGIDs(mortar_m, lag_mult_dof_map_);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> aux_mortar_m_trans =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(solidblock->RowMap(), 81, false));
  aux_mortar_m_trans->Add(*aux_mortar_m, true, -1.0 * (1.0 - solid_time_int_param), 0.0);
  aux_mortar_m_trans->Complete(*lag_mult_dof_map_, solidblock->RangeMap(), true);

  mat.Assign(0, 3, Core::LinAlg::View, *aux_mortar_m_trans);

  // ---------Addressing contribution to block (6,4)
  Teuchos::RCP<Core::LinAlg::SparseMatrix> aux_mortar_d =
      Mortar::MatrixRowTransformGIDs(mortar_d, lag_mult_dof_map_);

  aux_mortar_d->Scale(-1.0 / fluid_timescale);
  aux_mortar_d->Complete(fluidblock->FullDomainMap(), *lag_mult_dof_map_, true);

  mat.Assign(3, 1, Core::LinAlg::View, *aux_mortar_d);

  // ---------Addressing contribution to block (4,6)
  aux_mortar_d = Mortar::MatrixRowTransformGIDs(mortar_d, lag_mult_dof_map_);
  Teuchos::RCP<Core::LinAlg::SparseMatrix> aux_mortar_d_trans =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(fluidblock->FullRowMap(), 81, false));
  aux_mortar_d_trans->Add(
      *aux_mortar_d, true, 1.0 * (1.0 - fluid_time_int_param) / fluid_res_scale, 0.0);
  aux_mortar_d_trans->Complete(*lag_mult_dof_map_, fluidblock->FullRangeMap(), true);

  mat.Assign(1, 3, Core::LinAlg::View, *aux_mortar_d_trans);

  /*--------------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block
  if (fluid_shape_deriv != Teuchos::null)
  {
    const Core::Adapter::Coupling& coup_fluid_ale = fluid_ale_coupling();

    // extract submatrices
    Core::LinAlg::SparseMatrix& fluid_mesh_inner_inner = fluid_shape_deriv->Matrix(0, 0);
    Core::LinAlg::SparseMatrix& fluid_mesh_interf_inner = fluid_shape_deriv->Matrix(1, 0);
    Core::LinAlg::SparseMatrix& fluid_mesh_inner_interf = fluid_shape_deriv->Matrix(0, 1);
    Core::LinAlg::SparseMatrix& fluid_mesh_interf_interf = fluid_shape_deriv->Matrix(1, 1);

    // Adressing contribution to block (3,4)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> aux_fluid_mesh_inner_interf =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(fluid_mesh_inner_interf.RowMap(), 81, false));
    aux_fluid_mesh_inner_interf->Add(fluid_mesh_inner_interf, false, 1.0, 0.0);
    aux_fluid_mesh_inner_interf->Complete(
        fluidblock->DomainMap(), aux_fluid_mesh_inner_interf->RangeMap(), true);
    aux_fluidblock->Add(*aux_fluid_mesh_inner_interf, false, 1. / fluid_timescale, 1.0);

    // Adressing contribution to block (3,5)
    (*fluid_mesh_inner_inner_transform_)(fluid_shape_deriv->FullRowMap(),
        fluid_shape_deriv->FullColMap(), fluid_mesh_inner_inner, 1.,
        Core::Adapter::CouplingSlaveConverter(coup_fluid_ale), mat.Matrix(1, 2), false);

    // Adressing contribution to block (4,4)
    Teuchos::RCP<Core::LinAlg::SparseMatrix> aux_fluid_mesh_interf_interf =
        Teuchos::rcp(new Core::LinAlg::SparseMatrix(fluid_mesh_interf_interf.RowMap(), 81, false));
    aux_fluid_mesh_interf_interf->Add(fluid_mesh_interf_interf, false, 1.0, 0.0);
    aux_fluid_mesh_interf_interf->Complete(
        fluidblock->DomainMap(), aux_fluid_mesh_interf_interf->RangeMap(), true);
    aux_fluidblock->Add(*aux_fluid_mesh_interf_interf, false, 1. / fluid_timescale, 1.0);

    // Adressing contribution to block (4,5)
    (*fluid_mesh_inner_inner_transform_)(fluid_shape_deriv->FullRowMap(),
        fluid_shape_deriv->FullColMap(), fluid_mesh_interf_inner, 1.,
        Core::Adapter::CouplingMasterConverter(coup_fluid_ale), mat.Matrix(1, 2), false);
  }

  // finally assign fluid matrix to block (1,1)
  mat.Assign(1, 1, Core::LinAlg::View, *aux_fluidblock);

  // done. make sure all blocks are filled.
  mat.Complete();

  // Finally, take care of Dirichlet boundary conditions
  mat.ApplyDirichlet(*(dbcmaps_->CondMap()), true);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::scale_system(
    Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm =
      static_cast<bool>(Core::UTILS::IntegralValue<int>(fsimono, "INFNORMSCALING"));

  if (scaling_infnorm)
  {
    // do scaling of structure rows
    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);
    if (A->LeftScale(*srowsum_) or A->RightScale(*scolsum_) or
        mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 3).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(3, 0).EpetraMatrix()->RightScale(*scolsum_))
      FOUR_C_THROW("structure scaling failed");

    // do scaling of ale rows
    A = mat.Matrix(2, 2).EpetraMatrix();
    arowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    acolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*arowsum_);
    A->InvColSums(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2, 1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2, 3).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1, 2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(3, 2).EpetraMatrix()->RightScale(*acolsum_))
      FOUR_C_THROW("ale scaling failed");

    // do scaling of structure and ale rhs vectors
    Teuchos::RCP<Epetra_Vector> sx = extractor().ExtractVector(b, 0);
    Teuchos::RCP<Epetra_Vector> ax = extractor().ExtractVector(b, 2);

    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ax->Multiply(1.0, *arowsum_, *ax, 0.0)) FOUR_C_THROW("ale scaling failed");

    extractor().InsertVector(*sx, 0, b);
    extractor().InsertVector(*ax, 2, b);
  }
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::unscale_solution(
    Core::LinAlg::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm =
      static_cast<bool>(Core::UTILS::IntegralValue<int>(fsimono, "INFNORMSCALING"));

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_Vector> sy = extractor().ExtractVector(x, 0);
    Teuchos::RCP<Epetra_Vector> ay = extractor().ExtractVector(x, 2);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ay->Multiply(1.0, *acolsum_, *ay, 0.0)) FOUR_C_THROW("ale scaling failed");

    extractor().InsertVector(*sy, 0, x);
    extractor().InsertVector(*ay, 2, x);

    Teuchos::RCP<Epetra_Vector> sx = extractor().ExtractVector(b, 0);
    Teuchos::RCP<Epetra_Vector> ax = extractor().ExtractVector(b, 2);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0)) FOUR_C_THROW("structure scaling failed");
    if (ax->ReciprocalMultiply(1.0, *arowsum_, *ax, 0.0)) FOUR_C_THROW("ale scaling failed");

    extractor().InsertVector(*sx, 0, b);
    extractor().InsertVector(*ax, 2, b);

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);
    if (A->LeftScale(*srowsum_) or A->RightScale(*scolsum_) or
        mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 3).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(3, 0).EpetraMatrix()->RightScale(*scolsum_))
      FOUR_C_THROW("structure scaling failed");

    A = mat.Matrix(2, 2).EpetraMatrix();
    arowsum_->Reciprocal(*arowsum_);
    acolsum_->Reciprocal(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2, 1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2, 3).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1, 2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(3, 2).EpetraMatrix()->RightScale(*acolsum_))
      FOUR_C_THROW("ale scaling failed");
  }

  Epetra_Vector r(b.Map());
  mat.Apply(x, r);
  r.Update(1., b, 1.);

  Teuchos::RCP<Epetra_Vector> sr = extractor().ExtractVector(r, 0);
  Teuchos::RCP<Epetra_Vector> fr = extractor().ExtractVector(r, 1);
  Teuchos::RCP<Epetra_Vector> ar = extractor().ExtractVector(r, 2);
  Teuchos::RCP<Epetra_Vector> lmr = extractor().ExtractVector(r, 3);

  // increment additional ale residual
  aleresidual_->Update(-1., *ar, 0.);

  std::ios_base::fmtflags flags = utils()->out().flags();

  double n, ns, nf, na, nlm;
  r.Norm2(&n);
  sr->Norm2(&ns);
  fr->Norm2(&nf);
  ar->Norm2(&na);
  lmr->Norm2(&nlm);
  utils()->out() << std::scientific << "\nlinear solver quality:\n"
                 << "L_2-norms:\n"
                 << "   |r|=" << n << "   |rs|=" << ns << "   |rf|=" << nf << "   |ra|=" << na
                 << "   |rlm|=" << nlm << "\n";
  r.NormInf(&n);
  sr->NormInf(&ns);
  fr->NormInf(&nf);
  ar->NormInf(&na);
  lmr->NormInf(&nlm);
  utils()->out() << "L_inf-norms:\n"
                 << "   |r|=" << n << "   |rs|=" << ns << "   |rf|=" << nf << "   |ra|=" << na
                 << "   |rlm|=" << nlm << "\n";

  utils()->out().flags(flags);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::Evaluate(
    Teuchos::RCP<const Epetra_Vector> step_increment)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MortarMonolithicFluidSplitSaddlePoint::Evaluate");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // check whether all fields have the same time step size
  check_if_dts_same();
#endif

  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  Teuchos::RCP<const Epetra_Vector> ax;
  Teuchos::RCP<const Epetra_Vector> lagx;

  if (step_increment != Teuchos::null)
  {
    extract_field_vectors(step_increment, sx, fx, ax, lagx);
  }

  // Call all elements and assemble rhs and matrices
  // We only need the rhs here because NOX will ask for the rhs
  // only. But the Jacobian is stored internally and will be returned
  // later on without looking at x again!

  if (verbosity_ >= Inpar::FSI::verbosity_medium) utils()->out() << "\nEvaluate elements\n";

  {
    Teuchos::Time ts("structure", true);
    structure_field()->Evaluate(sx);
    if (verbosity_ >= Inpar::FSI::verbosity_medium)
      utils()->out() << "structure           : " << ts.totalElapsedTime(true) << " sec\n";
  }

  {
    Teuchos::Time ta("ale", true);
    ale_field()->Evaluate(ax);
    if (verbosity_ >= Inpar::FSI::verbosity_medium)
      utils()->out() << "ale                 : " << ta.totalElapsedTime(true) << " sec\n";
  }

  // transfer the current ale mesh positions to the fluid field
  Teuchos::RCP<Epetra_Vector> fluiddisp = ale_to_fluid(ale_field()->Dispnp());
  fluid_field()->apply_mesh_displacement(fluiddisp);

  {
    Teuchos::Time tf("fluid", true);
    fluid_field()->Evaluate(fx);
    if (verbosity_ >= Inpar::FSI::verbosity_medium)
      utils()->out() << "fluid                : " << tf.totalElapsedTime(true) << " sec\n";
  }

  {
    if (lagx != Teuchos::null)
    {
      Teuchos::Time tlm("lag_mult", true);
      lag_mult_->Update(1.0, *lag_mult_old_, 1.0, *lagx, 0.0);
      if (verbosity_ >= Inpar::FSI::verbosity_medium)
        utils()->out() << "Lagrange multiplier: " << tlm.totalElapsedTime(true) << " sec\n";
    }
  }

  if (verbosity_ >= Inpar::FSI::verbosity_medium) utils()->out() << "\n";
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::extract_field_vectors(
    Teuchos::RCP<const Epetra_Vector> x, Teuchos::RCP<const Epetra_Vector>& sx,
    Teuchos::RCP<const Epetra_Vector>& fx, Teuchos::RCP<const Epetra_Vector>& ax,
    Teuchos::RCP<const Epetra_Vector>& lagx)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MortarMonolithicFluidSplitSaddlePoint::extract_field_vectors");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  if (ddgpred_ == Teuchos::null)
    FOUR_C_THROW("Vector 'ddgpred_' has not been initialized properly.");
#endif

  // get the mortar structure to fluid coupling matrix M
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortar_m =
      coupling_solid_fluid_mortar_->GetMortarMatrixM();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mortar_m_transf =
      Mortar::MatrixRowTransformGIDs(mortar_m, lag_mult_dof_map_);

  // get the mortar fluid to structure coupling matrix D
  const Teuchos::RCP<Core::LinAlg::SparseMatrix> mortar_d =
      coupling_solid_fluid_mortar_->GetMortarMatrixD();
  Teuchos::RCP<Core::LinAlg::SparseMatrix> mortar_d_transf =
      Mortar::MatrixRowTransformGIDs(mortar_d, lag_mult_dof_map_);

  // ---------------------------------------------------------------------------
  // process structure unknowns
  // ---------------------------------------------------------------------------
  sx = extractor().ExtractVector(x, 0);
  // extract structure solution increment from NOX increment

  // ---------------------------------------------------------------------------
  // process fluid unknowns
  // ---------------------------------------------------------------------------
  // extract fluid solution increment from NOX increment
  fx = extractor().ExtractVector(x, 1);

  // ---------------------------------------------------------------------------
  // process ale unknowns
  // ---------------------------------------------------------------------------
  // extract inner ALE solution increment from NOX increment
  Teuchos::RCP<const Epetra_Vector> aox = extractor().ExtractVector(x, 2);

  // convert fluid interface velocities into ALE interface displacements
  Teuchos::RCP<Epetra_Vector> fcx = fluid_field()->Interface()->ExtractFSICondVector(fx);
  fluid_field()->velocity_to_displacement(fcx);
  Teuchos::RCP<Epetra_Vector> acx = fluid_to_ale_interface(fcx);

  // put inner and interface ALE solution increments together
  Teuchos::RCP<Epetra_Vector> a = ale_field()->Interface()->InsertOtherVector(aox);
  ale_field()->Interface()->InsertFSICondVector(acx, a);

  lagx = extractor().ExtractVector(x, 3);

  ax = a;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::update()
{
  // save Lagrange multiplier for the next time step
  lag_mult_old_->Update(1.0, *lag_mult_, 0.0);

  // call Update()-routine in base class to handle the single fields
  FSI::BlockMonolithic::update();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::output()
{
  structure_field()->Output();
  fluid_field()->Output();

  // output Lagrange multiplier
  OutputLambda();

  ale_field()->Output();

  if (structure_field()->get_constraint_manager()->HaveMonitor())
  {
    structure_field()->get_constraint_manager()->compute_monitor_values(
        structure_field()->Dispnp());
    if (comm_.MyPID() == 0) structure_field()->get_constraint_manager()->PrintMonitorValues();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::OutputLambda()
{
  /* 'lag_mult_' is only defined on the interface. So, insert 'lag_mult_' into
   * 'lambdafull' that is defined on the entire fluid field. Then, we need to write
   * output or restart data.
   */
  auto copy = Teuchos::rcp(new Epetra_Vector(*lag_mult_));
  copy->ReplaceMap(*fluid_field()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Vector> lambdafull = fluid_field()->Interface()->InsertFSICondVector(copy);
  const int uprestart = timeparams_.get<int>("RESTARTEVRY");
  const int upres = timeparams_.get<int>("RESULTSEVRY");
  if ((uprestart != 0 and fluid_field()->Step() % uprestart == 0) or
      (upres != 0 and fluid_field()->Step() % upres == 0))
    fluid_field()->DiscWriter()->write_vector("fsilambda", lambdafull);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::read_restart(int step)
{
  structure_field()->read_restart(step);
  fluid_field()->read_restart(step);

  // read Lagrange multiplier into fluid map
  Teuchos::RCP<Epetra_Vector> lambdafull =
      Teuchos::rcp(new Epetra_Vector(*fluid_field()->dof_row_map(), true));
  Core::IO::DiscretizationReader reader = Core::IO::DiscretizationReader(
      fluid_field()->discretization(), Global::Problem::Instance()->InputControlFile(), step);
  reader.read_vector(lambdafull, "fsilambda");
  auto lag_mult_old_on_fluid_map = fluid_field()->Interface()->ExtractFSICondVector(lambdafull);

  // Convert Lagrange multipliers to their actual map
  lag_mult_old_on_fluid_map->ReplaceMap(*lag_mult_dof_map_);
  lag_mult_old_ = Teuchos::RCP(new Epetra_Vector(*lag_mult_old_on_fluid_map));

  // Note: the above is normally enough. However, we can use the restart in order to periodically
  // repeat the fsi simulation (see AC-FS3I)
  lag_mult_ = Teuchos::RCP(new Epetra_Vector(*lag_mult_old_on_fluid_map));

  SetupSystem();

  ale_field()->read_restart(step);

  SetTimeStep(fluid_field()->Time(), fluid_field()->Step());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double FSI::MortarMonolithicFluidSplitSaddlePoint::select_dt_error_based() const
{
  // get time step size suggestions
  const double dtstr = get_ada_str_dt();           // based on all structure DOFs
  const double dtstrfsi = get_ada_str_fsi_dt();    // based on structure FSI DOFs
  const double dtflinner = get_ada_fl_inner_dt();  // based on inner fluid DOFs

  double dt = Dt();

  // select time step size based on error estimation
  if (is_ada_structure() and is_ada_fluid())
    dt = std::min(std::min(dtstr, dtstrfsi), dtflinner);
  else if (is_ada_structure() and (not is_ada_fluid()))
    dt = std::min(dtstr, dtstrfsi);
  else if ((not is_ada_structure()) and is_ada_fluid())
    dt = dtflinner;
  else
  {
    // no change in time step size based on structure or fluid field error estimation
  }

  return dt;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool FSI::MortarMonolithicFluidSplitSaddlePoint::set_accepted() const
{
  // get error norms
  const double strnorm = get_ada_strnorm();            // based on all structure DOFs
  const double strfsinorm = get_ada_str_fs_inorm();    // based on structure FSI DOFs
  const double flinnernorm = get_ada_fl_inner_norm();  // based on inner fluid DOFs

  bool accepted = std::max(strnorm, strfsinorm) < errtolstr_ and flinnernorm < errtolfl_;

  // in case error estimation in the fluid field is turned off:
  if (not is_ada_fluid()) accepted = std::max(strnorm, strfsinorm) < errtolstr_;

  // in case error estimation in the structure field is turned off:
  if (not is_ada_structure()) accepted = flinnernorm < errtolfl_;

  // no error based time adaptivity
  if ((not is_ada_structure()) and (not is_ada_fluid())) accepted = true;

  return accepted;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::create_node_owner_relationship(
    std::map<int, int>* nodeOwner, std::map<int, std::list<int>>* inverseNodeOwner,
    std::map<int, Core::Nodes::Node*>* fluidnodesPtr,
    std::map<int, Core::Nodes::Node*>* structuregnodesPtr,
    Teuchos::RCP<Core::FE::Discretization> structuredis,
    Teuchos::RCP<Core::FE::Discretization> fluiddis, const Inpar::FSI::Redistribute domain)
{
  FOUR_C_THROW("Not implemented, yet.");
}

FOUR_C_NAMESPACE_CLOSE
