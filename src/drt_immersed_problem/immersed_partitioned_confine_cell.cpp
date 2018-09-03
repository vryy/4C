/*!----------------------------------------------------------------------
\file immersed_partitioned_confine_cell.cpp

\brief partitioned immersed cell-ecm interaction via cell confinement in ecm pore

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

\level 3
*----------------------------------------------------------------------*/
#include "immersed_partitioned_confine_cell.H"

#include "../linalg/linalg_utils.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_poroelast/poro_scatra_base.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_adapter/ad_str_multiphysicswrapper_cellmigration.H"


IMMERSED::ImmersedPartitionedConfineCell::ImmersedPartitionedConfineCell(
    const Teuchos::ParameterList& params, const Epetra_Comm& comm)
    : ImmersedPartitioned(comm)
{
  // get pointer to fluid search tree from ParameterList
  fluid_SearchTree_ = params.get<Teuchos::RCP<GEO::SearchTree>>("RCPToFluidSearchTree");

  // get pointer to cell search tree from ParameterList
  cell_SearchTree_ = params.get<Teuchos::RCP<GEO::SearchTree>>("RCPToCellSearchTree");

  // get pointer to the current position map of the cell
  currpositions_cell_ =
      params.get<std::map<int, LINALG::Matrix<3, 1>>*>("PointerToCurrentPositionsCell");

  // get pointer to the current position map of the cell
  currpositions_ECM_ =
      params.get<std::map<int, LINALG::Matrix<3, 1>>*>("PointerToCurrentPositionsECM");

  // get pointer to cell structure
  Teuchos::RCP<ADAPTER::MultiphysicsStructureWrapperCellMigration> multiphysicswrapper =
      params.get<Teuchos::RCP<ADAPTER::MultiphysicsStructureWrapperCellMigration>>(
          "RCPToCellStructure");

  if (multiphysicswrapper == Teuchos::null)
    dserror("no pointer to MultiphysicsStructureWrapperCellMigration provided");

  cellstructure_ = multiphysicswrapper->GetFSIStructureWrapperPtr();

  // create instance of poroelast subproblem
  poroscatra_subproblem_ = params.get<Teuchos::RCP<POROELAST::PoroScatraBase>>("RCPToPoroScatra");

  // create instance of poroelast subproblem
  poroscatra_subproblem_ = params.get<Teuchos::RCP<POROELAST::PoroScatraBase>>("RCPToPoroScatra");

  // set pointer to poro fpsi structure
  porostructure_ = poroscatra_subproblem_->PoroField()->StructureField();

  // check object pointers
  if (fluid_SearchTree_ == Teuchos::null) dserror("no pointer to fluid_SearchTree_ provided !");
  if (cell_SearchTree_ == Teuchos::null) dserror("no pointer to cell_SearchTree_ provided !");
  if (currpositions_cell_ == NULL) dserror("no pointer to currpositions_cell_ provided !");
  if (currpositions_ECM_ == NULL) dserror("no pointer to currpositions_ECM_ provided !");
  if (cellstructure_ == Teuchos::null) dserror("no pointer to cellstructure_ provided !");
  if (poroscatra_subproblem_ == Teuchos::null)
    dserror("no pointer to poroscatra_subproblem_ provided !");

  // important variables for parallel simulations
  myrank_ = comm.MyPID();
  numproc_ = comm.NumProc();

  // get pointer to global problem
  globalproblem_ = DRT::Problem::Instance();

  backgroundfluiddis_ = globalproblem_->GetDis("porofluid");
  backgroundstructuredis_ = globalproblem_->GetDis("structure");
  immerseddis_ = globalproblem_->GetDis("cell");

  // construct immersed exchange manager. singleton class that makes immersed variables comfortably
  // accessible from everywhere in the code
  exchange_manager_ = DRT::ImmersedFieldExchangeManager::Instance();
  exchange_manager_->SetIsPureConfinementSimulation(
      params.get<bool>("IsPureConfinementSimulation") == true);
  exchange_manager_->InitializePorosityAtGPMap();
  exchange_manager_->SetIsInitialized(true);

  // get coupling variable
  displacementcoupling_ = globalproblem_->CellMigrationParams()
                              .sublist("CONFINEMENT MODULE")
                              .get<std::string>("COUPVARIABLE") == "Displacement";
  if (displacementcoupling_ and myrank_ == 0)
    std::cout
        << "\n Coupling variable for partitioned Cell-ECM Confinement scheme :  Displacements "
        << std::endl;
  else if (!displacementcoupling_ and myrank_ == 0)
    std::cout << "\n Coupling variable for partitioned Cell-ECM Confinement scheme :  Force "
              << std::endl;

  // vector of penalty traction on cell bdry int points integrated over cell surface
  cell_penalty_traction_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()), true));
  // force increment
  cell_delta_penalty_traction_ =
      Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()), true));
  // gap integrated over surface
  penalty_gap_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()), true));
  // current nodal normals
  curr_nodal_normals_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()), true));

  exchange_manager_->SetPointerToECMPenaltyTraction(cell_penalty_traction_);
  exchange_manager_->SetPointerToGap(penalty_gap_);
  exchange_manager_->SetPointerToCurrentNodalNormals(curr_nodal_normals_);

  // get pointer to cell search tree from ParameterList
  cell_SearchTree_ = params.get<Teuchos::RCP<GEO::SearchTree>>("RCPToCellSearchTree");

  // PSEUDO2D switch
  isPseudo2D_ = DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(), "PSEUDO2D");

  // penalty constraint enforcement parameters
  penalty_start_ = globalproblem_->CellMigrationParams()
                       .sublist("CONFINEMENT MODULE")
                       .get<double>("PENALTY_START");
  penalty_init_ = globalproblem_->CellMigrationParams()
                      .sublist("CONFINEMENT MODULE")
                      .get<double>("PENALTY_INIT");

  // initialize the parameter initialize_cell_ which determines whether or not first time step is
  // pre-simulation
  initialize_cell_ =
      DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(), "INITIALIZE_CELL");
  timestep_ = Dt();
  if (initialize_cell_)
    initial_timestep_ = globalproblem_->CellMigrationParams().get<double>("INITIAL_TIMESTEP");
  else
    initial_timestep_ = timestep_;

  // get number of timesteps of cell initialization
  initialization_steps_ = globalproblem_->CellMigrationParams().get<int>("INITIALIZATION_STEPS");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int IMMERSED::ImmersedPartitionedConfineCell::Init(const Teuchos::ParameterList& params)
{
  // reset the setup flag
  SetIsSetup(false);

  // do all init stuff here

  // set isinit_ flag true
  SetIsInit(true);

  return 0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedConfineCell::Setup()
{
  // make sure Init(...) was called first
  CheckIsInit();

  // get parameters for nox
  const Teuchos::ParameterList& noxparams =
      globalproblem_->CellMigrationParams().sublist("CONFINEMENT MODULE");
  SetDefaultParameters(noxparams, NOXParameterList());
  // noxparameterlist_.print();

  // set flag issetup true
  SetIsSetup(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedConfineCell::CouplingOp(
    const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  ReinitTransferVectors();

  ResetImmersedInformation();

  // DISPLACEMENT COUPLING
  if (displacementcoupling_)
  {
    const Teuchos::RCP<Epetra_Vector> idispn = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL ImmersedOp
    ////////////////////
    PrepareImmersedOp();
    Teuchos::RCP<Epetra_Vector> idispnp = ImmersedOp(
        cellstructure_->Interface()->ExtractIMMERSEDCondVector(cell_penalty_traction_), fillFlag);

    ////////////////////
    // CALL BackgroundOp
    ////////////////////
    PrepareBackgroundOp();
    BackgroundOp(Teuchos::null, fillFlag);

    int err = F.Update(1.0, *idispnp, -1.0, *idispn, 0.0);
    if (err != 0) dserror("Vector update of Coupling-residual returned err=%d", err);
  }
  // FORCE COUPLING
  else if (!displacementcoupling_)
  {
    dserror("force coupling not tested, yet");
    const Teuchos::RCP<Epetra_Vector> cell_reactforcen = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL BackgroundOp
    ////////////////////
    PrepareBackgroundOp();
    BackgroundOp(cell_reactforcen, fillFlag);

    ////////////////////
    // CALL ImmersedOp
    ////////////////////


    int err = F.Update(1.0, *cell_penalty_traction_, -1.0, *cell_reactforcen, 0.0);

    if (err != 0) dserror("Vector update of FSI-residual returned err=%d", err);

  }  // displacement / force coupling

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedConfineCell::BackgroundOp(
    Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values, const FillType fillFlag)
{
  IMMERSED::ImmersedPartitioned::BackgroundOp(backgrd_dirichlet_values, fillFlag);

  if (fillFlag == User)
  {
    dserror("fillFlag == User : not yet implemented");
  }
  else
  {
    if (myrank_ == 0)
      std::cout << "BackgroundOp is empty. So far, this is only a one way coupled Problem.\n"
                   "It is assumed, that the cell does not deform the ECM when it is compressed.\n"
                << std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedConfineCell::ImmersedOp(
    Teuchos::RCP<Epetra_Vector> bdry_traction, const FillType fillFlag)
{
  IMMERSED::ImmersedPartitioned::ImmersedOp(bdry_traction, fillFlag);

  if (fillFlag == User)
  {
    dserror("fillFlag == User : not yet implemented");
    return Teuchos::null;
  }
  else
  {
    // prescribe neumann values at structural boundary dofs
    cellstructure_->ApplyImmersedInterfaceForces(Teuchos::null, bdry_traction);

    // solve cell
    cellstructure_->Solve();

    return cellstructure_->ExtractImmersedInterfaceDispnp();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedConfineCell::PrepareBackgroundOp()
{
  // do nothing so far

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedConfineCell::PrepareImmersedOp()
{
  immerseddis_->SetState(0, "displacement", cellstructure_->Dispnp());
  backgroundfluiddis_->SetState(0, "dispnp", poroscatra_subproblem_->FluidField()->Dispnp());
  backgroundstructuredis_->SetState(0, "displacement", porostructure_->Dispnp());

  double structsearchradiusfac =
      DRT::Problem::Instance()->ImmersedMethodParams().get<double>("STRCT_SRCHRADIUS_FAC");

  // determine subset of fluid discretization which is potentially underlying the immersed
  // discretization
  //
  // get state
  Teuchos::RCP<const Epetra_Vector> celldisplacements = cellstructure_->Dispnp();

  // find current positions for immersed structural discretization
  std::map<int, LINALG::Matrix<3, 1>> my_currpositions_cell;
  for (int lid = 0; lid < immerseddis_->NumMyRowNodes(); ++lid)
  {
    const DRT::Node* node = immerseddis_->lRowNode(lid);
    LINALG::Matrix<3, 1> currpos;
    std::vector<int> dofstoextract(3);
    std::vector<double> mydisp(3);

    // get the current displacement
    immerseddis_->Dof(node, 0, dofstoextract);
    DRT::UTILS::ExtractMyValues(*celldisplacements, mydisp, dofstoextract);

    currpos(0) = node->X()[0] + mydisp.at(0);
    currpos(1) = node->X()[1] + mydisp.at(1);
    currpos(2) = node->X()[2] + mydisp.at(2);

    my_currpositions_cell[node->Id()] = currpos;
  }

  // communicate local currpositions:
  // map with current cell positions should be same on all procs
  // to make use of the advantages of ghosting the cell redundantly
  // on all procs.
  int procs[Comm().NumProc()];
  for (int i = 0; i < Comm().NumProc(); i++) procs[i] = i;
  LINALG::Gather<int, LINALG::Matrix<3, 1>>(
      my_currpositions_cell, *currpositions_cell_, Comm().NumProc(), &procs[0], Comm());

  // get bounding box of current configuration of structural dis
  const LINALG::Matrix<3, 2> structBox = GEO::getXAABBofDis(*immerseddis_, *currpositions_cell_);
  double max_radius =
      sqrt(pow(structBox(0, 0) - structBox(0, 1), 2) + pow(structBox(1, 0) - structBox(1, 1), 2) +
           pow(structBox(2, 0) - structBox(2, 1), 2));
  // search for background elements within a certain radius around the center of the immersed
  // bounding box
  LINALG::Matrix<3, 1> boundingboxcenter;
  boundingboxcenter(0) = structBox(0, 0) + (structBox(0, 1) - structBox(0, 0)) * 0.5;
  boundingboxcenter(1) = structBox(1, 0) + (structBox(1, 1) - structBox(1, 0)) * 0.5;
  boundingboxcenter(2) = structBox(2, 0) + (structBox(2, 1) - structBox(2, 0)) * 0.5;

  // search background elements covered by bounding box of cell
  curr_subset_of_backgrounddis_ = fluid_SearchTree_->searchElementsInRadius(*backgroundfluiddis_,
      *currpositions_ECM_, boundingboxcenter, structsearchradiusfac * max_radius, 0);

  if (curr_subset_of_backgrounddis_.empty() == false)
    std::cout << "\nPrepareBackgroundOp returns "
              << curr_subset_of_backgrounddis_.begin()->second.size()
              << " background elements on Proc " << Comm().MyPID() << std::endl;

  if (myrank_ == 0) std::cout << "   Update Immersed Information ..." << std::endl;

  Teuchos::ParameterList params;

  if (curr_subset_of_backgrounddis_.empty() == false)
    EvaluateImmersedNoAssembly(params, backgroundfluiddis_, &curr_subset_of_backgrounddis_,
        cell_SearchTree_, currpositions_cell_, (int)FLD::update_immersed_information);
  else
  {
    // do nothing : a proc without subset of backgrd dis. does not need to enter evaluation,
    //              since cell is ghosted on all procs and no communication has to be performed.
  }

  // calculate compressive force from ECM to cell
  CalcECMTractionOnStructure();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedConfineCell::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  if (myrank_ == 0) std::cout << "Cell Predictor: " << std::endl;
  cellstructure_->PrepareTimeStep();
  if (myrank_ == 0) std::cout << "Poro Predictor: " << std::endl;
  poroscatra_subproblem_->PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedConfineCell::InitialGuess()
{
  if (displacementcoupling_)
    return cellstructure_->PredictImmersedInterfaceDispnp();
  else  // FORCE COUPLING
  {
    return cellstructure_->Interface()->ExtractIMMERSEDCondVector(cell_penalty_traction_);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedConfineCell::CalcECMTractionOnStructure()
{
  double penalty = penalty_start_;

  Teuchos::ParameterList params;
  params.set<std::string>("action", "calc_ecm_traction");
  params.set<std::string>("backgrddisname", "structure");
  params.set<std::string>("immerseddisname", "cell");

  double normofstructbdrytraction = -1234.0;
  double normofdeltastructbdrytraction = -1234.0;

  DRT::AssembleStrategy cell_str_bdry_strategy(0,  // struct dofset for row
      0,                                           // struct dofset for column
      Teuchos::null,                               // matrix 1
      Teuchos::null,                               //
      cell_delta_penalty_traction_,                // vector 1
      penalty_gap_,                                // vector 2
      Teuchos::null                                //
  );
  if (myrank_ == 0)
  {
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
    std::cout << "###   Calculate force from ECM transmitted to cell surface due to confinement    "
                 "               "
              << std::endl;
  }

  exchange_manager_->SetGapMax(-1234.0);
  exchange_manager_->SetVoidMax(-1234.0);
  exchange_manager_->SetGapMin(1234.0);
  exchange_manager_->SetVoidMin(1234.0);
  exchange_manager_->SetDeltaPorosityMax(0.0);

  if (!initialize_cell_)
    penalty = penalty_start_;  //*sqrt(((IterationCounter())[0]));
  else
    penalty = penalty_init_;

  params.set<double>("penalty", penalty);
  params.set<int>("during_init", initialize_cell_);
  EvaluateInterpolationCondition(
      immerseddis_, params, cell_str_bdry_strategy, "IMMERSEDCoupling", -1);

  // Apply relaxed ecm force increment
  cell_penalty_traction_->Update(1.0, *cell_delta_penalty_traction_,
      0.0);  // cell_delta_penalty_traction contains full force, thus this*0.0
  cell_delta_penalty_traction_->Norm2(&normofdeltastructbdrytraction);

  normofstructbdrytraction = 0.0;
  cell_penalty_traction_->Norm2(&normofstructbdrytraction);
  if (myrank_ == 0)
  {
    std::cout << "###   Penalty Parameter:  " << penalty << std::endl;
    std::cout << "###   Max delta_phi:  " << std::setprecision(6)
              << exchange_manager_->GetDeltaPorosityMax() << std::endl;
    std::cout << "###   MAX void = " << std::setprecision(6) << exchange_manager_->GetVoidMax()
              << "    MIN void = " << std::setprecision(6) << exchange_manager_->GetVoidMin()
              << std::endl;
    std::cout << "###   MAX gap  = " << std::setprecision(6) << exchange_manager_->GetGapMax()
              << "   MIN gap  = " << std::setprecision(6) << exchange_manager_->GetGapMin()
              << std::endl;
    std::cout << "###   Max gap origin:    [" << std::setprecision(6)
              << (*(exchange_manager_->GetMaxGapSpacePoint()))(0) << " "
              << (*(exchange_manager_->GetMaxGapSpacePoint()))(1) << " "
              << (*(exchange_manager_->GetMaxGapSpacePoint()))(2) << "]" << std::endl;
    std::cout << "###   Max gap end point: ["
              << (*(exchange_manager_->GetMaxGapSearchDirection()))(0) << " "
              << (*(exchange_manager_->GetMaxGapSearchDirection()))(1) << " "
              << (*(exchange_manager_->GetMaxGapSearchDirection()))(2) << "]" << std::endl;
    std::cout << "###   Min gap origin:    [" << std::setprecision(6)
              << (*(exchange_manager_->GetMinGapSpacePoint()))(0) << " "
              << (*(exchange_manager_->GetMinGapSpacePoint()))(1) << " "
              << (*(exchange_manager_->GetMinGapSpacePoint()))(2) << "]" << std::endl;
    std::cout << "###   Min gap end point: ["
              << (*(exchange_manager_->GetMinGapSearchDirection()))(0) << " "
              << (*(exchange_manager_->GetMinGapSearchDirection()))(1) << " "
              << (*(exchange_manager_->GetMinGapSearchDirection()))(2) << "]" << std::endl;
    std::cout << "###   Norm of Boundary ECM interaction force vector:    " << std::setprecision(13)
              << normofstructbdrytraction << std::endl;
    std::cout << "###   Norm of Boundary ECM interaction force increment: " << std::setprecision(13)
              << normofdeltastructbdrytraction << std::endl;
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedConfineCell::ReadRestart(int step)
{
  cellstructure_->ReadRestart(step);
  poroscatra_subproblem_->ReadRestart(step);
  SetTimeStep(poroscatra_subproblem_->PoroField()->Time(), step);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedConfineCell::PrepareOutput()
{
  poroscatra_subproblem_->PrepareOutput();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedConfineCell::Output()
{
  cellstructure_->Output();
  poroscatra_subproblem_->Output();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedConfineCell::Update()
{
  cellstructure_->Update();
  poroscatra_subproblem_->Update();
  exchange_manager_->UpdatePorosityAtGPOldTimestep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedConfineCell::ResetImmersedInformation()
{
  Teuchos::ParameterList params;
  params.set<int>("action", FLD::reset_immersed_ele);
  params.set<int>("Physical Type", poroscatra_subproblem_->FluidField()->PhysicalType());
  backgroundfluiddis_->Evaluate(params);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedConfineCell::SetFieldDt()
{
  if (initialize_cell_ and Step() == 0 and (IterationCounter())[0] == 0)
  {
    if (myrank_ == 0)
    {
      std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
      std::cout << "++++++                 CELL INITIALIZATION STEP                   ++++++"
                << std::endl;
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
      std::cout << "dt=" << initial_timestep_ << "\n" << std::endl;
    }
    SetDt(initial_timestep_);
    poroscatra_subproblem_->SetDt(initial_timestep_);
    poroscatra_subproblem_->PoroField()->SetDt(initial_timestep_);
    poroscatra_subproblem_->FluidField()->SetDt(initial_timestep_);
    porostructure_->SetDt(initial_timestep_);
    poroscatra_subproblem_->ScaTraField()->SetDt(initial_timestep_);
    cellstructure_->SetDt(initial_timestep_);
  }
  else if (initialize_cell_ and Step() == initialization_steps_)
  {
    if (myrank_ == 0)
    {
      std::cout << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
      std::cout << "++++++              CELL INITIALIZATION SUCCESSFUL                ++++++"
                << std::endl;
      std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << std::endl;
      std::cout << " RESET TIMESTEP TO dt=" << timestep_ << "\n" << std::endl;
    }
    SetDt(timestep_);
    poroscatra_subproblem_->SetDt(timestep_);
    poroscatra_subproblem_->PoroField()->SetDt(timestep_);
    poroscatra_subproblem_->FluidField()->SetDt(timestep_);
    porostructure_->SetDt(timestep_);
    poroscatra_subproblem_->ScaTraField()->SetDt(timestep_);
    cellstructure_->SetDt(timestep_);
    initialize_cell_ = 0;

    double simpleECMInteractionConstant =
        ((((exchange_manager_->GetPointerToPenaltyTractionAtGPMap()->find(3))->second).at(0))
                .Norm2()) /
        0.2;
    exchange_manager_->SetSimpleECMInteractionConstant(simpleECMInteractionConstant);
    std::cout << "Set Slope of ECM Interaction Force to " << simpleECMInteractionConstant
              << std::endl;
  }

  return;
}
