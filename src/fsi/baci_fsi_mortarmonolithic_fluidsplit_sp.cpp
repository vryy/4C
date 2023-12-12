/*----------------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problem with non-matching grids using a monolithic scheme
in saddle-point formulation with Lagrange multipliers discretized on the fluid interface

\level 1
*/

/*----------------------------------------------------------------------------*/


#include "baci_fsi_mortarmonolithic_fluidsplit_sp.H"

#include "baci_adapter_ale_fsi.H"
#include "baci_adapter_fld_fluid_fsi.H"
#include "baci_adapter_str_fsiwrapper.H"
#include "baci_ale_utils_mapextractor.H"
#include "baci_constraint_manager.H"
#include "baci_coupling_adapter.H"
#include "baci_coupling_adapter_converter.H"
#include "baci_coupling_adapter_mortar.H"
#include "baci_fluid_utils_mapextractor.H"
#include "baci_fsi_debugwriter.H"
#include "baci_fsi_statustest.H"
#include "baci_io.H"
#include "baci_io_control.H"
#include "baci_io_pstream.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_blocksparsematrix.H"
#include "baci_linalg_mapextractor.H"
#include "baci_linalg_matrixtransform.H"
#include "baci_linalg_sparsematrix.H"
#include "baci_linalg_utils_sparse_algebra_print.H"
#include "baci_mortar_utils.H"
#include "baci_structure_aux.H"
#include "baci_utils_exceptions.H"

#include <Epetra_Comm.h>
#include <Teuchos_TimeMonitor.hpp>

BACI_NAMESPACE_OPEN

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
  intersectionmaps.push_back(FluidField()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(FluidField()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmap =
      CORE::LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

  // Check whether the intersection is empty
  if (intersectionmap->NumGlobalElements() != 0)
  {
    //      std::cout << "Slave interface nodes with Dirichlet boundary condition "
    //                "(input file numbering):" << std::endl;
    //      for (int i=0; i < (int)FluidField()->Discretization()->NumMyRowNodes(); i++)
    //      {
    //        // get all nodes and add them
    //        int gid = FluidField()->Discretization()->NodeRowMap()->GID(i);

    //        // do only nodes that I have in my discretization
    //        if (!FluidField()->Discretization()->NodeColMap()->MyGID(gid)) continue;
    //        DRT::Node* node = FluidField()->Discretization()->gNode(gid);
    //        if (!node) dserror("Cannot find node with gid %",gid);

    //        std::vector<int> nodedofs = FluidField()->Discretization()->Dof(node);

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

    dserror(errormsg.str());
  }
  // ---------------------------------------------------------------------------

  notsetup_ = true;

  coupling_solid_fluid_mortar_ = Teuchos::rcp(new CORE::ADAPTER::CouplingMortar(
      DRT::Problem::Instance()->NDim(), DRT::Problem::Instance()->MortarCouplingParams(),
      DRT::Problem::Instance()->ContactDynamicParams(),
      DRT::Problem::Instance()->SpatialApproximationType()));

  ale_inner_interf_transform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fluid_mesh_inner_inner_transform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);

  CreateLagrangeMultiplierDofRowMap();
  SetLagMult();

#ifdef DEBUG
  if (coupling_solid_fluid_mortar_ == Teuchos::null)
  {
    dserror("Allocation of 'coupling_solid_fluid_mortar_' failed.");
  }
  if (ale_inner_interf_transform_ == Teuchos::null)
  {
    dserror("Allocation of 'ale_inner_interf_transform_' failed.");
  }
  if (fluid_mesh_inner_inner_transform_ == Teuchos::null)
  {
    dserror("Allocation of 'fluid_mesh_inner_inner_transform_' failed.");
  }
  if (lag_mult_ == Teuchos::null)
  {
    dserror("Allocation of 'lag_mult_' failed.");
  }
#endif
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::SetLagMult()
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
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

    SetDefaultParameters(fsidyn, NOXParameterList());

    // we use non-matching meshes at the interface
    // mortar with: structure = master, fluid = slave

    const int ndim = DRT::Problem::Instance()->NDim();

    // get coupling objects
    CORE::ADAPTER::Coupling& interface_coup_fluid_ale = InterfaceFluidAleCoupling();

    /* structure to fluid
     * coupling condition at the fsi interface:
     * displacements (=number spacial dimensions) are coupled
     * e.g.: 3D: coupleddof = [1, 1, 1]
     */
    std::vector<int> coupleddof(ndim, 1);

    coupling_solid_fluid_mortar_->Setup(StructureField()->Discretization(),
        FluidField()->Discretization(), AleField()->WriteAccessDiscretization(), coupleddof,
        "FSICoupling", comm_, true);

    // fluid to ale at the interface
    interface_coup_fluid_ale.SetupConditionCoupling(*FluidField()->Discretization(),
        FluidField()->Interface()->FSICondMap(), *AleField()->Discretization(),
        AleField()->Interface()->FSICondMap(), "FSICoupling", ndim);

    CORE::ADAPTER::Coupling& coup_fluid_ale = FluidAleCoupling();

    // the fluid-ale coupling always matches
    const Epetra_Map* fluidnodemap = FluidField()->Discretization()->NodeRowMap();
    const Epetra_Map* alenodemap = AleField()->Discretization()->NodeRowMap();

    coup_fluid_ale.SetupCoupling(*FluidField()->Discretization(), *AleField()->Discretization(),
        *fluidnodemap, *alenodemap, ndim);

    FluidField()->SetMeshMap(coup_fluid_ale.MasterDofMap());

    CreateCombinedDofRowMap();

    /*------------------------------------------------------------------------*/
    // Switch fluid to interface split block matrix
    FluidField()->UseBlockMatrix(true);

    // build ale system matrix in splitted system
    AleField()->CreateSystemMatrix(AleField()->Interface());

    aleresidual_ = Teuchos::rcp(new Epetra_Vector(*AleField()->Interface()->OtherMap()));

    // -------------------------------------------------------------------------
    // Build the global Dirichlet map extractor
    SetupDBCMapExtractor();
    // -------------------------------------------------------------------------#

    // enable debugging
    if (DRT::INPUT::IntegralValue<int>(fsidyn, "DEBUGOUTPUT") & 2)
    {
      pcdbg_ = Teuchos::rcp(new UTILS::MonolithicDebugWriter(*this));
    }

    CreateSystemMatrix();

    notsetup_ = false;
  }

  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    const bool restartfrompartfsi =
        DRT::INPUT::IntegralValue<bool>(timeparams_, "RESTART_FROM_PART_FSI");
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
void FSI::MortarMonolithicFluidSplitSaddlePoint::CreateLagrangeMultiplierDofRowMap()
{
  const int num_glob_elem_fluid_interface =
      FluidField()->Interface()->FSICondMap()->NumGlobalElements();
  const int num_loc_elem_fluid_interface = FluidField()->Interface()->FSICondMap()->NumMyElements();
  const int max_gid_ale = AleField()->DofRowMap()->MaxAllGID();
  lag_mult_dof_map_ = Teuchos::rcp(new Epetra_Map(
      num_glob_elem_fluid_interface, num_loc_elem_fluid_interface, max_gid_ale + 1, comm_));
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::CreateCombinedDofRowMap()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(StructureField()->DofRowMap());
  vecSpaces.push_back(FluidField()->DofRowMap());
  vecSpaces.push_back(AleField()->Interface()->OtherMap());
  vecSpaces.push_back(lag_mult_dof_map_);

  if (vecSpaces[1]->NumGlobalElements() == 0)
    dserror("No inner fluid equations. Splitting not possible.");

  SetDofRowMaps(vecSpaces);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::CreateSystemMatrix()
{
  FSI::BlockMonolithic::CreateSystemMatrix(systemmatrix_, false);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Combo> FSI::MortarMonolithicFluidSplitSaddlePoint::CreateStatusTest(
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
      "DISPL residual", Extractor(), 0, nlParams.get<double>("Tol dis res L2"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "DISPL residual", Extractor(), 0, nlParams.get<double>("Tol dis res Inf"),
      ::NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("DISPL update", Extractor(), 0,
          nlParams.get<double>("Tol dis inc L2"), ::NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("DISPL update", Extractor(), 0,
          nlParams.get<double>("Tol dis inc Inf"), ::NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  AddStatusTest(structureDisp_L2);

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
  fluidvel.push_back(FluidField()->InnerVelocityRowMap());
  fluidvel.push_back(Teuchos::null);
  CORE::LINALG::MultiMapExtractor fluidvelextract(*DofRowMap(), fluidvel);

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
  AddStatusTest(innerFluidVel_L2);

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
  fluidpress.push_back(FluidField()->PressureRowMap());
  fluidpress.push_back(Teuchos::null);
  CORE::LINALG::MultiMapExtractor fluidpressextract(*DofRowMap(), fluidpress);

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
  AddStatusTest(fluidPress_L2);

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
      "LAGMULT residual", Extractor(), 3, nlParams.get<double>("Tol fsi res L2"),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormF> lag_mult_inf = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "LAGMULT residual", Extractor(), 3, nlParams.get<double>("Tol fsi res Inf"),
      ::NOX::Abstract::Vector::MaxNorm, NOX::FSI::PartialNormF::Unscaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> lag_mult_update_L2 =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("LAGMULT update", Extractor(), 3,
          nlParams.get<double>("Tol fsi inc L2"), ::NOX::Abstract::Vector::TwoNorm,
          NOX::FSI::PartialNormUpdate::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> lag_mult_update_inf =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("LAGMULT update", Extractor(), 3,
          nlParams.get<double>("Tol fsi inc Inf"), ::NOX::Abstract::Vector::MaxNorm,
          NOX::FSI::PartialNormUpdate::Unscaled));

  // tests needed to adapt relative tolerance of the linear solver
  AddStatusTest(lag_mult_L2);

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
void FSI::MortarMonolithicFluidSplitSaddlePoint::SetupDBCMapExtractor()
{
  /* Dirichlet maps for structure and fluid do not intersect with interface map.
   * ALE Dirichlet map might intersect with interface map, but ALE interface
   * DOFs are not part of the final system of equations. Hence, we just need the
   * intersection of inner ALE DOFs with Dirichlet ALE DOFs.
   */
  std::vector<Teuchos::RCP<const Epetra_Map>> aleintersectionmaps;
  aleintersectionmaps.push_back(AleField()->GetDBCMapExtractor()->CondMap());
  aleintersectionmaps.push_back(AleField()->Interface()->OtherMap());
  Teuchos::RCP<const Epetra_Map> aleintersectionmap =
      CORE::LINALG::MultiMapExtractor::IntersectMaps(aleintersectionmaps);

  // Merge Dirichlet maps of structure, fluid and ALE to global FSI Dirichlet map
  std::vector<Teuchos::RCP<const Epetra_Map>> dbcmaps;
  dbcmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
  dbcmaps.push_back(FluidField()->GetDBCMapExtractor()->CondMap());
  dbcmaps.push_back(aleintersectionmap);

  Teuchos::RCP<const Epetra_Map> dbcmap = CORE::LINALG::MultiMapExtractor::MergeMaps(dbcmaps);

  // Finally, create the global FSI Dirichlet map extractor
  dbcmaps_ = Teuchos::rcp(new CORE::LINALG::MapExtractor(*DofRowMap(), dbcmap, true));
  if (dbcmaps_ == Teuchos::null) dserror("Creation of FSI Dirichlet map extractor failed.");
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>
FSI::MortarMonolithicFluidSplitSaddlePoint::SystemMatrix() const
{
  return systemmatrix_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::InitialGuess(
    Teuchos::RCP<Epetra_Vector> initial_guess)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MortarMonolithicFluidSplitSaddlePoint::InitialGuess");

  Teuchos::RCP<const Epetra_Vector> lag_mult_initial_guess =
      Teuchos::rcp(new Epetra_Vector(*lag_mult_dof_map_, true));

  CombineFieldVectors(*initial_guess, StructureField()->InitialGuess(),
      FluidField()->InitialGuess(), AleField()->InitialGuess(), lag_mult_initial_guess, true);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::CombineFieldVectors(Epetra_Vector& f,
    Teuchos::RCP<const Epetra_Vector> solid_vector, Teuchos::RCP<const Epetra_Vector> fluid_vector,
    Teuchos::RCP<const Epetra_Vector> ale_vector, Teuchos::RCP<const Epetra_Vector> lag_mult_vector,
    bool fullvectors)
{
  if (fullvectors)
  {
    // extract inner DOFs from slave vectors
    Teuchos::RCP<const Epetra_Vector> ale_other_vector =
        AleField()->Interface()->ExtractOtherVector(ale_vector);

    Extractor().AddVector(*solid_vector, 0, f);
    Extractor().AddVector(*fluid_vector, 1, f);
    Extractor().AddVector(*ale_other_vector, 2, f);
    Extractor().AddVector(*lag_mult_vector, 3, f);
  }
  else
  {
    Extractor().AddVector(*solid_vector, 0, f);
    Extractor().AddVector(*fluid_vector, 1, f);
    Extractor().AddVector(*ale_vector, 2, f);
    Extractor().AddVector(*lag_mult_vector, 3, f);
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::SetupRHSResidual(Epetra_Vector& f)
{
  // get single field residuals
  Teuchos::RCP<const Epetra_Vector> solid_single_field_rhs_vector =
      Teuchos::rcp(new Epetra_Vector(*StructureField()->RHS()));
  Teuchos::RCP<const Epetra_Vector> fluid_single_field_rhs_vector =
      Teuchos::rcp(new Epetra_Vector(*FluidField()->RHS()));
  Teuchos::RCP<const Epetra_Vector> ale_single_field_rhs_vector =
      Teuchos::rcp(new Epetra_Vector(*AleField()->RHS()));
  Teuchos::RCP<Epetra_Vector> lag_mult_rhs_vector =
      Teuchos::rcp(new Epetra_Vector(*lag_mult_dof_map_, true));

  // put the single field residuals together
  CombineFieldVectors(f, solid_single_field_rhs_vector, fluid_single_field_rhs_vector,
      ale_single_field_rhs_vector, lag_mult_rhs_vector, true);

  // add additional ale residual to avoid incremental ale errors
  Extractor().AddVector(*aleresidual_, 2, f);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::SetupRHSLambda(Epetra_Vector& f)
{
  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double solid_time_int_param = StructureField()->TimIntParam();
  const double fluid_time_int_param = FluidField()->TimIntParam();
  const double fluid_res_scale = FluidField()->ResidualScaling();

  // get the mortar structure to fluid coupling matrix M
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> mortar_m =
      coupling_solid_fluid_mortar_->GetMortarMatrixM();
  Teuchos::RCP<CORE::LINALG::SparseMatrix> mortar_m_transf =
      MORTAR::MatrixRowTransformGIDs(mortar_m, lag_mult_dof_map_);

  // get the mortar fluid to structure coupling matrix D
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> mortar_d =
      coupling_solid_fluid_mortar_->GetMortarMatrixD();
  Teuchos::RCP<CORE::LINALG::SparseMatrix> mortar_d_transf =
      MORTAR::MatrixRowTransformGIDs(mortar_d, lag_mult_dof_map_);

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
      StructureField()->Interface()->InsertFSICondVector(lag_mult_old_rhs_struc_interf);
  Teuchos::RCP<Epetra_Vector> lag_mult_old_rhs_fluid_interf_full =
      FluidField()->Interface()->InsertFSICondVector(lag_mult_old_rhs_fluid_interf);

  lag_mult_old_rhs_fluid_interf_full->Scale(-1.0 / fluid_res_scale);

  // add lagrange multiplier
  Extractor().AddVector(*lag_mult_old_rhs_struc_interf_full, 0, f);
  Extractor().AddVector(*lag_mult_old_rhs_fluid_interf_full, 1, f);

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
      StructureField()->Interface()->InsertFSICondVector(lag_mult_step_increment_rhs_struc_interf);
  Teuchos::RCP<Epetra_Vector> lag_mult_step_increment_rhs_fluid_interf_full =
      FluidField()->Interface()->InsertFSICondVector(lag_mult_step_increment_rhs_fluid_interf);

  lag_mult_step_increment_rhs_struc_interf_full->Scale(1.0 * (1. - solid_time_int_param));
  lag_mult_step_increment_rhs_fluid_interf_full->Scale(
      -1.0 * (1. - fluid_time_int_param) / fluid_res_scale);

  // add lagrange multiplier
  Extractor().AddVector(*lag_mult_step_increment_rhs_struc_interf_full, 0, f);
  Extractor().AddVector(*lag_mult_step_increment_rhs_fluid_interf_full, 1, f);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::SetupRHSFirstiter(Epetra_Vector& f)
{
  // old interface velocity of fluid field
  const Teuchos::RCP<const Epetra_Vector> fluid_veln = FluidField()->ExtractInterfaceVeln();

  // get the mortar structure to fluid coupling matrix M
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> mortar_m =
      coupling_solid_fluid_mortar_->GetMortarMatrixM();
  Teuchos::RCP<CORE::LINALG::SparseMatrix> mortar_m_transf =
      MORTAR::MatrixRowTransformGIDs(mortar_m, lag_mult_dof_map_);

  // get the mortar fluid to structure coupling matrix D
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> mortar_d =
      coupling_solid_fluid_mortar_->GetMortarMatrixD();
  Teuchos::RCP<CORE::LINALG::SparseMatrix> mortar_d_transf =
      MORTAR::MatrixRowTransformGIDs(mortar_d, lag_mult_dof_map_);

  // get fluid shape derivatives matrix
  const Teuchos::RCP<const CORE::LINALG::BlockSparseMatrixBase> fluid_shape_deriv =
      FluidField()->ShapeDerivatives();

  // get ale matrix
  const Teuchos::RCP<const CORE::LINALG::BlockSparseMatrixBase> aleblock =
      AleField()->BlockSystemMatrix();

  // extract ale submatrix
  const CORE::LINALG::SparseMatrix& ale_inner_interf = aleblock->Matrix(0, 1);

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
    const CORE::LINALG::SparseMatrix& fluid_mesh_inner_interf = fluid_shape_deriv->Matrix(0, 1);

    rhs = Teuchos::rcp(new Epetra_Vector(fluid_mesh_inner_interf.RangeMap(), true));

    fluid_mesh_inner_interf.Apply(*fluid_veln, *rhs);

    rhs->Scale(Dt());

    rhs = FluidField()->Interface()->InsertOtherVector(rhs);

    Extractor().AddVector(*rhs, 1, f);
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
    const CORE::LINALG::SparseMatrix& fluid_mesh_interf_interf = fluid_shape_deriv->Matrix(1, 1);
    rhs = Teuchos::rcp(new Epetra_Vector(fluid_mesh_interf_interf.RangeMap(), true));

    fluid_mesh_interf_interf.Apply(*fluid_veln, *rhs);
    rhs->Scale(Dt());
    rhs = FluidField()->Interface()->InsertFSICondVector(rhs);

    Extractor().AddVector(*rhs, 1, f);
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

  Extractor().AddVector(*rhs, 2, f);
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

  Extractor().AddVector(*rhs, 3, f);
  // ----------end of term 1

  // ----------addressing term 2
  rhs = Teuchos::rcp(new Epetra_Vector(*lag_mult_dof_map_, true));

  mortar_m_transf->Apply(*ddgpred_, *rhs);
  rhs->Scale(-1.);

  Extractor().AddVector(*rhs, 3, f);
  // ----------end of term 2
  // ----------end of lagrange multiplier
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::SetupSystemMatrix(
    CORE::LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MortarMonolithicFluidSplitSaddlePoint::SetupSystemMatrix");

  // get the mortar structure to fluid coupling matrix M
  const Teuchos::RCP<const CORE::LINALG::SparseMatrix> mortar_m =
      coupling_solid_fluid_mortar_->GetMortarMatrixM();

  // get the mortar fluid to structure coupling matrix D
  const Teuchos::RCP<const CORE::LINALG::SparseMatrix> mortar_d =
      coupling_solid_fluid_mortar_->GetMortarMatrixD();

  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double solid_time_int_param = StructureField()->TimIntParam();
  const double fluid_time_int_param = FluidField()->TimIntParam();
  const double fluid_res_scale = FluidField()->ResidualScaling();

  // time scaling factor for fluid
  const double fluid_timescale = FluidField()->TimeScaling();

  // get fluid shape derivatives matrix
  const Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> fluid_shape_deriv =
      FluidField()->ShapeDerivatives();

  // get single field block matrices
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> solidblock = StructureField()->SystemMatrix();
  const Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> fluidblock =
      FluidField()->BlockSystemMatrix();
  const Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> aleblock =
      AleField()->BlockSystemMatrix();

  // extract submatrices
  const CORE::LINALG::SparseMatrix& fluid_inner_inner = fluidblock->Matrix(0, 0);
  const CORE::LINALG::SparseMatrix& fluid_interf_inner = fluidblock->Matrix(1, 0);
  const CORE::LINALG::SparseMatrix& fluid_inner_interf = fluidblock->Matrix(0, 1);
  const CORE::LINALG::SparseMatrix& fluid_interf_interf = fluidblock->Matrix(1, 1);
  const CORE::LINALG::SparseMatrix& ale_inner_inner = aleblock->Matrix(0, 0);
  const CORE::LINALG::SparseMatrix& ale_inner_interf = aleblock->Matrix(0, 1);

  // ---------------------------------------------------------------------------
  // BEGIN building the global 6x6 system matrix
  // ---------------------------------------------------------------------------
  // Contributions to blocks in system matrix are listed separately.
  // Block numbering in comments ranges from (1,1) to (6,6).

  // ---------Addressing contribution to blocks (1,1),(1,2),(2,1),(2,2)
  mat.Assign(0, 0, CORE::LINALG::View, *solidblock);

  // ---------Addressing contribution to blocks (3,3),(3,4),(4,3),(4,4)
  Teuchos::RCP<CORE::LINALG::SparseMatrix> aux_fluidblock =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(fluidblock->FullRowMap(), 108, false));
  aux_fluidblock->Add(fluid_inner_inner, false, 1.0, 0.0);
  aux_fluidblock->Add(fluid_interf_inner, false, 1.0, 1.0);
  aux_fluidblock->Add(fluid_inner_interf, false, 1.0, 1.0);
  aux_fluidblock->Add(fluid_interf_interf, false, 1.0, 1.0);
  aux_fluidblock->Complete(fluidblock->FullDomainMap(), fluidblock->FullRangeMap(), true);

  // ---------Addressing contribution to block (5,4)
  Teuchos::RCP<CORE::LINALG::SparseMatrix> aux_ale_inner_interf =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(ale_inner_inner.RowMap(), 81, false));
  (*ale_inner_interf_transform_)(aleblock->FullRowMap(), aleblock->FullColMap(), ale_inner_interf,
      1., CORE::ADAPTER::CouplingSlaveConverter(InterfaceFluidAleCoupling()),
      *aux_ale_inner_interf);

  aux_ale_inner_interf->Scale(1. / fluid_timescale);
  aux_ale_inner_interf->Complete(fluidblock->DomainMap(), aux_ale_inner_interf->RangeMap(), true);

  mat.Assign(2, 1, CORE::LINALG::View, *aux_ale_inner_interf);

  // ---------Addressing contribution to block (5,5)
  mat.Assign(2, 2, CORE::LINALG::View, ale_inner_inner);

  // ---------Addressing contribution to block (6,2)
  Teuchos::RCP<CORE::LINALG::SparseMatrix> aux_mortar_m =
      MORTAR::MatrixRowTransformGIDs(mortar_m, lag_mult_dof_map_);
  aux_mortar_m->Complete(solidblock->DomainMap(), *lag_mult_dof_map_, true);

  mat.Assign(3, 0, CORE::LINALG::View, *aux_mortar_m);

  // ---------Addressing contribution to block (2,6)
  aux_mortar_m = MORTAR::MatrixRowTransformGIDs(mortar_m, lag_mult_dof_map_);
  Teuchos::RCP<CORE::LINALG::SparseMatrix> aux_mortar_m_trans =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(solidblock->RowMap(), 81, false));
  aux_mortar_m_trans->Add(*aux_mortar_m, true, -1.0 * (1.0 - solid_time_int_param), 0.0);
  aux_mortar_m_trans->Complete(*lag_mult_dof_map_, solidblock->RangeMap(), true);

  mat.Assign(0, 3, CORE::LINALG::View, *aux_mortar_m_trans);

  // ---------Addressing contribution to block (6,4)
  Teuchos::RCP<CORE::LINALG::SparseMatrix> aux_mortar_d =
      MORTAR::MatrixRowTransformGIDs(mortar_d, lag_mult_dof_map_);

  aux_mortar_d->Scale(-1.0 / fluid_timescale);
  aux_mortar_d->Complete(fluidblock->FullDomainMap(), *lag_mult_dof_map_, true);

  mat.Assign(3, 1, CORE::LINALG::View, *aux_mortar_d);

  // ---------Addressing contribution to block (4,6)
  aux_mortar_d = MORTAR::MatrixRowTransformGIDs(mortar_d, lag_mult_dof_map_);
  Teuchos::RCP<CORE::LINALG::SparseMatrix> aux_mortar_d_trans =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(fluidblock->FullRowMap(), 81, false));
  aux_mortar_d_trans->Add(
      *aux_mortar_d, true, 1.0 * (1.0 - fluid_time_int_param) / fluid_res_scale, 0.0);
  aux_mortar_d_trans->Complete(*lag_mult_dof_map_, fluidblock->FullRangeMap(), true);

  mat.Assign(1, 3, CORE::LINALG::View, *aux_mortar_d_trans);

  /*--------------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block
  if (fluid_shape_deriv != Teuchos::null)
  {
    const CORE::ADAPTER::Coupling& coup_fluid_ale = FluidAleCoupling();

    // extract submatrices
    CORE::LINALG::SparseMatrix& fluid_mesh_inner_inner = fluid_shape_deriv->Matrix(0, 0);
    CORE::LINALG::SparseMatrix& fluid_mesh_interf_inner = fluid_shape_deriv->Matrix(1, 0);
    CORE::LINALG::SparseMatrix& fluid_mesh_inner_interf = fluid_shape_deriv->Matrix(0, 1);
    CORE::LINALG::SparseMatrix& fluid_mesh_interf_interf = fluid_shape_deriv->Matrix(1, 1);

    // Adressing contribution to block (3,4)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> aux_fluid_mesh_inner_interf =
        Teuchos::rcp(new CORE::LINALG::SparseMatrix(fluid_mesh_inner_interf.RowMap(), 81, false));
    aux_fluid_mesh_inner_interf->Add(fluid_mesh_inner_interf, false, 1.0, 0.0);
    aux_fluid_mesh_inner_interf->Complete(
        fluidblock->DomainMap(), aux_fluid_mesh_inner_interf->RangeMap(), true);
    aux_fluidblock->Add(*aux_fluid_mesh_inner_interf, false, 1. / fluid_timescale, 1.0);

    // Adressing contribution to block (3,5)
    (*fluid_mesh_inner_inner_transform_)(fluid_shape_deriv->FullRowMap(),
        fluid_shape_deriv->FullColMap(), fluid_mesh_inner_inner, 1.,
        CORE::ADAPTER::CouplingSlaveConverter(coup_fluid_ale), mat.Matrix(1, 2), false);

    // Adressing contribution to block (4,4)
    Teuchos::RCP<CORE::LINALG::SparseMatrix> aux_fluid_mesh_interf_interf =
        Teuchos::rcp(new CORE::LINALG::SparseMatrix(fluid_mesh_interf_interf.RowMap(), 81, false));
    aux_fluid_mesh_interf_interf->Add(fluid_mesh_interf_interf, false, 1.0, 0.0);
    aux_fluid_mesh_interf_interf->Complete(
        fluidblock->DomainMap(), aux_fluid_mesh_interf_interf->RangeMap(), true);
    aux_fluidblock->Add(*aux_fluid_mesh_interf_interf, false, 1. / fluid_timescale, 1.0);

    // Adressing contribution to block (4,5)
    (*fluid_mesh_inner_inner_transform_)(fluid_shape_deriv->FullRowMap(),
        fluid_shape_deriv->FullColMap(), fluid_mesh_interf_inner, 1.,
        CORE::ADAPTER::CouplingMasterConverter(coup_fluid_ale), mat.Matrix(1, 2), false);
  }

  // finally assign fluid matrix to block (1,1)
  mat.Assign(1, 1, CORE::LINALG::View, *aux_fluidblock);

  // done. make sure all blocks are filled.
  mat.Complete();

  // Finally, take care of Dirichlet boundary conditions
  mat.ApplyDirichlet(*(dbcmaps_->CondMap()), true);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::Evaluate(
    Teuchos::RCP<const Epetra_Vector> step_increment)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MortarMonolithicFluidSplitSaddlePoint::Evaluate");

#ifdef DEBUG
  // check whether all fields have the same time step size
  CheckIfDtsSame();
#endif

  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  Teuchos::RCP<const Epetra_Vector> ax;
  Teuchos::RCP<const Epetra_Vector> lagx;

  if (step_increment != Teuchos::null)
  {
    ExtractFieldVectors(step_increment, sx, fx, ax, lagx);
  }

  // Call all elements and assemble rhs and matrices
  // We only need the rhs here because NOX will ask for the rhs
  // only. But the Jacobian is stored internally and will be returned
  // later on without looking at x again!

  if (verbosity_ >= INPAR::FSI::verbosity_medium) Utils()->out() << "\nEvaluate elements\n";

  {
    Teuchos::Time ts("structure", true);
    StructureField()->Evaluate(sx);
    if (verbosity_ >= INPAR::FSI::verbosity_medium)
      Utils()->out() << "structure           : " << ts.totalElapsedTime(true) << " sec\n";
  }

  {
    Teuchos::Time ta("ale", true);
    AleField()->Evaluate(ax);
    if (verbosity_ >= INPAR::FSI::verbosity_medium)
      Utils()->out() << "ale                 : " << ta.totalElapsedTime(true) << " sec\n";
  }

  // transfer the current ale mesh positions to the fluid field
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluid(AleField()->Dispnp());
  FluidField()->ApplyMeshDisplacement(fluiddisp);

  {
    Teuchos::Time tf("fluid", true);
    FluidField()->Evaluate(fx);
    if (verbosity_ >= INPAR::FSI::verbosity_medium)
      Utils()->out() << "fluid                : " << tf.totalElapsedTime(true) << " sec\n";
  }

  {
    if (lagx != Teuchos::null)
    {
      Teuchos::Time tlm("lag_mult", true);
      lag_mult_->Update(1.0, *lag_mult_old_, 1.0, *lagx, 0.0);
      if (verbosity_ >= INPAR::FSI::verbosity_medium)
        Utils()->out() << "Lagrange multiplier: " << tlm.totalElapsedTime(true) << " sec\n";
    }
  }

  if (verbosity_ >= INPAR::FSI::verbosity_medium) Utils()->out() << "\n";
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::ExtractFieldVectors(
    Teuchos::RCP<const Epetra_Vector> x, Teuchos::RCP<const Epetra_Vector>& sx,
    Teuchos::RCP<const Epetra_Vector>& fx, Teuchos::RCP<const Epetra_Vector>& ax,
    Teuchos::RCP<const Epetra_Vector>& lagx)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MortarMonolithicFluidSplitSaddlePoint::ExtractFieldVectors");

#ifdef DEBUG
  if (ddgpred_ == Teuchos::null) dserror("Vector 'ddgpred_' has not been initialized properly.");
#endif

  // get the mortar structure to fluid coupling matrix M
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> mortar_m =
      coupling_solid_fluid_mortar_->GetMortarMatrixM();
  Teuchos::RCP<CORE::LINALG::SparseMatrix> mortar_m_transf =
      MORTAR::MatrixRowTransformGIDs(mortar_m, lag_mult_dof_map_);

  // get the mortar fluid to structure coupling matrix D
  const Teuchos::RCP<CORE::LINALG::SparseMatrix> mortar_d =
      coupling_solid_fluid_mortar_->GetMortarMatrixD();
  Teuchos::RCP<CORE::LINALG::SparseMatrix> mortar_d_transf =
      MORTAR::MatrixRowTransformGIDs(mortar_d, lag_mult_dof_map_);

  // ---------------------------------------------------------------------------
  // process structure unknowns
  // ---------------------------------------------------------------------------
  sx = Extractor().ExtractVector(x, 0);
  // extract structure solution increment from NOX increment

  // ---------------------------------------------------------------------------
  // process fluid unknowns
  // ---------------------------------------------------------------------------
  // extract fluid solution increment from NOX increment
  fx = Extractor().ExtractVector(x, 1);

  // ---------------------------------------------------------------------------
  // process ale unknowns
  // ---------------------------------------------------------------------------
  // extract inner ALE solution increment from NOX increment
  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x, 2);

  // convert fluid interface velocities into ALE interface displacements
  Teuchos::RCP<Epetra_Vector> fcx = FluidField()->Interface()->ExtractFSICondVector(fx);
  FluidField()->VelocityToDisplacement(fcx);
  Teuchos::RCP<Epetra_Vector> acx = FluidToAleInterface(fcx);

  // put inner and interface ALE solution increments together
  Teuchos::RCP<Epetra_Vector> a = AleField()->Interface()->InsertOtherVector(aox);
  AleField()->Interface()->InsertFSICondVector(acx, a);

  lagx = Extractor().ExtractVector(x, 3);

  ax = a;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::Update()
{
  // save Lagrange multiplier for the next time step
  lag_mult_old_->Update(1.0, *lag_mult_, 0.0);

  // call Update()-routine in base class to handle the single fields
  FSI::BlockMonolithic::Update();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::Output()
{
  StructureField()->Output();
  FluidField()->Output();

  // output Lagrange multiplier
  OutputLambda();

  AleField()->Output();

  if (StructureField()->GetConstraintManager()->HaveMonitor())
  {
    StructureField()->GetConstraintManager()->ComputeMonitorValues(StructureField()->Dispnp());
    if (comm_.MyPID() == 0) StructureField()->GetConstraintManager()->PrintMonitorValues();
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
  copy->ReplaceMap(*FluidField()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Vector> lambdafull = FluidField()->Interface()->InsertFSICondVector(copy);
  const int uprestart = timeparams_.get<int>("RESTARTEVRY");
  const int upres = timeparams_.get<int>("RESULTSEVRY");
  if ((uprestart != 0 and FluidField()->Step() % uprestart == 0) or
      (upres != 0 and FluidField()->Step() % upres == 0))
    FluidField()->DiscWriter()->WriteVector("fsilambda", lambdafull);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::ReadRestart(int step)
{
  StructureField()->ReadRestart(step);
  FluidField()->ReadRestart(step);

  // read Lagrange multiplier into fluid map
  Teuchos::RCP<Epetra_Vector> lambdafull =
      Teuchos::rcp(new Epetra_Vector(*FluidField()->DofRowMap(), true));
  IO::DiscretizationReader reader = IO::DiscretizationReader(FluidField()->Discretization(), step);
  reader.ReadVector(lambdafull, "fsilambda");
  auto lag_mult_old_on_fluid_map = FluidField()->Interface()->ExtractFSICondVector(lambdafull);

  // Convert Lagrange multipliers to their actual map
  lag_mult_old_on_fluid_map->ReplaceMap(*lag_mult_dof_map_);
  lag_mult_old_ = Teuchos::RCP(new Epetra_Vector(*lag_mult_old_on_fluid_map));

  // Note: the above is normally enough. However, we can use the restart in order to periodically
  // repeat the fsi simulation (see AC-FS3I)
  lag_mult_ = Teuchos::RCP(new Epetra_Vector(*lag_mult_old_on_fluid_map));

  SetupSystem();

  AleField()->ReadRestart(step);

  SetTimeStep(FluidField()->Time(), FluidField()->Step());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double FSI::MortarMonolithicFluidSplitSaddlePoint::SelectDtErrorBased() const
{
  // get time step size suggestions
  const double dtstr = GetAdaStrDt();          // based on all structure DOFs
  const double dtstrfsi = GetAdaStrFSIDt();    // based on structure FSI DOFs
  const double dtflinner = GetAdaFlInnerDt();  // based on inner fluid DOFs

  double dt = Dt();

  // select time step size based on error estimation
  if (IsAdaStructure() and IsAdaFluid())
    dt = std::min(std::min(dtstr, dtstrfsi), dtflinner);
  else if (IsAdaStructure() and (not IsAdaFluid()))
    dt = std::min(dtstr, dtstrfsi);
  else if ((not IsAdaStructure()) and IsAdaFluid())
    dt = dtflinner;
  else
  {
    // no change in time step size based on structure or fluid field error estimation
  }

  return dt;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool FSI::MortarMonolithicFluidSplitSaddlePoint::SetAccepted() const
{
  // get error norms
  const double strnorm = GetAdaStrnorm();          // based on all structure DOFs
  const double strfsinorm = GetAdaStrFSInorm();    // based on structure FSI DOFs
  const double flinnernorm = GetAdaFlInnerNorm();  // based on inner fluid DOFs

  bool accepted = std::max(strnorm, strfsinorm) < errtolstr_ and flinnernorm < errtolfl_;

  // in case error estimation in the fluid field is turned off:
  if (not IsAdaFluid()) accepted = std::max(strnorm, strfsinorm) < errtolstr_;

  // in case error estimation in the structure field is turned off:
  if (not IsAdaStructure()) accepted = flinnernorm < errtolfl_;

  // no error based time adaptivity
  if ((not IsAdaStructure()) and (not IsAdaFluid())) accepted = true;

  return accepted;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::MortarMonolithicFluidSplitSaddlePoint::CreateNodeOwnerRelationship(
    std::map<int, int>* nodeOwner, std::map<int, std::list<int>>* inverseNodeOwner,
    std::map<int, DRT::Node*>* fluidnodesPtr, std::map<int, DRT::Node*>* structuregnodesPtr,
    Teuchos::RCP<DRT::Discretization> structuredis, Teuchos::RCP<DRT::Discretization> fluiddis,
    const INPAR::FSI::Redistribute domain)
{
  dserror("Not implemented, yet.");
}

BACI_NAMESPACE_CLOSE
