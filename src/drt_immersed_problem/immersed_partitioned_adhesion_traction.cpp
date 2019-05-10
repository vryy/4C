/*!----------------------------------------------------------------------

\brief partitioned immersed cell-ecm interaction via adhesion traction

\level 2

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/
#include "immersed_partitioned_adhesion_traction.H"

#include "str_model_evaluator_partitioned_ssi_adhesiondynamics.H"

#include "../drt_structure/stru_aux.H"
#include "../drt_poroelast/poro_scatra_base.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_structure_new/str_dbc.H"
#include "../drt_structure_new/str_model_evaluator_multiphysics.H"
#include "../drt_adapter/ad_str_multiphysicswrapper_cellmigration.H"
#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../linalg/linalg_utils.H"
#include "../drt_ssi/ssi_partitioned_2wc.H"


IMMERSED::ImmersedPartitionedAdhesionTraction::ImmersedPartitionedAdhesionTraction(
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

  // get pointer to the current position map of the ECM
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

  // get pointer structure-scatra interaction (ssi) subproblem
  cellscatra_subproblem_ = params.get<Teuchos::RCP<SSI::SSI_Part2WC>>("RCPToCellScatra");

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
  if (cellscatra_subproblem_ == Teuchos::null)
    dserror("no pointer to cellscatra_subproblem_ provided !");

  // important variables for parallel simulations
  myrank_ = comm.MyPID();
  numproc_ = comm.NumProc();

  // get pointer to global problem
  globalproblem_ = DRT::Problem::Instance();

  // get pointers to discretizations
  backgroundfluiddis_ = globalproblem_->GetDis("porofluid");
  backgroundstructuredis_ = globalproblem_->GetDis("structure");
  immerseddis_ = globalproblem_->GetDis("cell");
  immersedscatradis_ = globalproblem_->GetDis("cellscatra");

  // construct immersed exchange manager. singleton class that makes immersed variables comfortably
  // accessible from everywhere in the code
  exchange_manager_ = DRT::ImmersedFieldExchangeManager::Instance();
  exchange_manager_->SetIsPureAdhesionSimulation(true);

  // get adhesion condition
  dirichletcoupling_ = globalproblem_->CellMigrationParams()
                           .sublist("ADHESION MODULE")
                           .get<std::string>("COUPMETHOD") == "Dirichlet";
  if (dirichletcoupling_ and myrank_ == 0)
    std::cout << " Coupling method for partitioned Cell-ECM Adhesion Dynamics scheme :  Dirichlet "
              << std::endl;
  else if (!dirichletcoupling_ and myrank_ == 0)
    std::cout << " Coupling method for partitioned Cell-ECM Adhesion Dynamics scheme :  Penalty "
              << std::endl;

  if (dirichletcoupling_)
  {
    // vector of adhesion forces in ecm
    cell_adhesion_disp_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()), true));
    Freact_cell_ = cellstructure_->Freact();
  }
  else  // violation coupling via penalty method
  {
    // get adhesion force pointer from model evaluator and set pointer in immersed exchange manager
    cell_adhesion_forces_ =
        Teuchos::rcp_dynamic_cast<::STR::MODELEVALUATOR::PartitionedSSIAdhesionDynamics>(
            multiphysicswrapper->CellmigrationModelEvaluator()->GetModelEvaluatorFromMap(
                STR::MODELEVALUATOR::MultiphysicType::mt_ssi))
            ->GetAdhesionForcesPtr();
    exchange_manager_->SetPointerCellAdhesionForce(cell_adhesion_forces_);

    // set pointer to adhesion fixpoint coordinates in model evaluator and exchange manager
    cell_adhesion_nod_coords_ = Teuchos::rcp(new Epetra_Vector(*(immerseddis_->DofRowMap()), true));
    Teuchos::rcp_dynamic_cast<::STR::MODELEVALUATOR::PartitionedSSIAdhesionDynamics>(
        multiphysicswrapper->CellmigrationModelEvaluator()->GetModelEvaluatorFromMap(
            STR::MODELEVALUATOR::MultiphysicType::mt_ssi))
        ->SetAdhesionFixpointPtr(cell_adhesion_nod_coords_);

    // get adhesion reset concentration
    reset_conc_ =
        globalproblem_->CellMigrationParams().sublist("ADHESION MODULE").get<double>("RESET_CONC");
    double scatra_tol = globalproblem_->CellMigrationParams()
                            .sublist("SCALAR TRANSPORT")
                            .sublist("NONLINEAR")
                            .get<double>("CONVTOL");
    if (scatra_tol > reset_conc_)
      dserror("Scatra tolerance higher than adhesion reset concentration!");
  }

  // set pointer to the adhesion force distributed to the ecm
  ecm_adhesion_forces_ = Teuchos::rcp(
      new Epetra_Vector(*(poroscatra_subproblem_->StructureField()->DofRowMap()), true));
  exchange_manager_->SetPointerECMAdhesionForce(ecm_adhesion_forces_);

  // PSEUDO2D switch
  isPseudo2D_ = DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(), "PSEUDO2D");

  // method for treatment of artificial ecm solid phase
  artificial_ecm_treatment_ = DRT::INPUT::IntegralValue<int>(
      globalproblem_->CellMigrationParams(), "ARTIFICIAL_ECM_METHOD");

  // initial invalid immersed information
  immersed_information_invalid_ = true;

  // output after every fixed-point iteration?
  output_evry_nlniter_ = (DRT::INPUT::IntegralValue<int>(
      globalproblem_->ImmersedMethodParams(), "OUTPUT_EVRY_NLNITER"));

  // switch deciding whether adhesion state has to be (re-)determined
  determine_adhesion_information_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int IMMERSED::ImmersedPartitionedAdhesionTraction::Init(const Teuchos::ParameterList& params)
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
void IMMERSED::ImmersedPartitionedAdhesionTraction::Setup()
{
  // make sure Init(...) was called first
  CheckIsInit();

  // get parameters for nox
  const Teuchos::ParameterList& noxparams =
      globalproblem_->CellMigrationParams().sublist("ADHESION MODULE");
  SetDefaultParameters(noxparams, NOXParameterList());
  // noxparameterlist_.print();

  // set flag issetup true
  SetIsSetup(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::CouplingOp(
    const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  ReinitTransferVectors();

  // reset immersed information, if penalty method is applied
  // this is not necessary when using fixed Dirichlet values as the adhering points do not change
  if (!dirichletcoupling_) ResetImmersedInformation();

  // get old reaction force
  const Teuchos::RCP<Epetra_Vector> cell_reactforcen = Teuchos::rcp(new Epetra_Vector(x));

  ////////////////////
  // CALL BackgroundOp
  ////////////////////
  PrepareBackgroundOp();
  BackgroundOp(cell_reactforcen, fillFlag);

  ////////////////////
  // CALL ImmersedOp
  ////////////////////
  ImmersedOp(Teuchos::null, fillFlag);

  int err = -1;

  if (dirichletcoupling_)
    err = F.Update(1.0, *Freact_cell_, -1.0, *cell_reactforcen, 0.0);
  else
    err = F.Update(1.0, *cell_adhesion_forces_, -1.0, *cell_reactforcen, 0.0);

  if (err != 0) dserror("Vector update of FSI-residual returned err=%d", err);

  // write output after every solve of ECM and Cell
  // current limitations:
  // max 100 partitioned iterations and max 100 timesteps in total
  if (output_evry_nlniter_)
  {
    int iter = ((IterationCounter())[0]);
    poroscatra_subproblem_->PrepareOutput();
    poroscatra_subproblem_->FluidField()->Output(
        (Step() * 100) + (iter - 1), Time() - Dt() * ((100 - iter) / 100.0));
    cellstructure_->PrepareOutput();
    cellstructure_->Output(
        false, (Step() * 100) + (iter - 1), Time() - Dt() * ((100 - iter) / 100.0));
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::BackgroundOp(
    Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values, const FillType fillFlag)
{
  IMMERSED::ImmersedPartitioned::BackgroundOp(backgrd_dirichlet_values, fillFlag);

  if (fillFlag == User)
  {
    dserror("fillFlag == User : not yet implemented");
  }
  else
  {
    Teuchos::ParameterList params;

    if (immersed_information_invalid_)
    {
      if (myrank_ == 0) std::cout << "   \nUpdate Immersed Information ...\n" << std::endl;

      if (curr_subset_of_backgrounddis_.empty() == false)
        EvaluateImmersedNoAssembly(params, backgroundfluiddis_, &curr_subset_of_backgrounddis_,
            cell_SearchTree_, currpositions_cell_, (int)FLD::update_immersed_information);
      else
      {
        // do nothing : a proc without subset of backgrd dis. does not need to enter evaluation,
        //              since cell is ghosted on all procs and no communication has to be performed.
      }

      // do only if we want to apply a method for the cancellation of the effect from
      // the artificial ecm solid phase onto the physical ecm solid phase.
      if (artificial_ecm_treatment_)
      {
        ///////////////////////////////////////////////////////////////////////////////////
        // prepare for removing the effect of the ecm mechanics in the artificial domain
        ///////////////////////////////////////////////////////////////////////////////////
        if (myrank_ == 0)
          std::cout << "\n   Determine ECM-Cell Pseudo-Boundary and label elements ...\n"
                    << std::endl;
        std::vector<int> ElesInFirstImmersedRow;
        std::vector<int> nodes =
            DetermineElesInFirstImmersedRow(backgroundfluiddis_, &ElesInFirstImmersedRow);
        CreateVolumeCondition(backgroundstructuredis_, nodes, DRT::Condition::CellFocalAdhesion,
            "ImmersedAdhesionElements", true);
        for (int id = 0; id < (int)ElesInFirstImmersedRow.size(); ++id)
          std::cout << "PROC " << myrank_ << ": " << ElesInFirstImmersedRow[id] << std::endl;
        Comm().Barrier();

        exchange_manager_->SetEles(ElesInFirstImmersedRow);
      }

      immersed_information_invalid_ = false;
    }  // if immersed_information_invalid_

    // distribute reaction force from cell boundary nodes to adjacent ECM nodes
    DistributeAdhesionForce(backgrd_dirichlet_values);
    immersed_information_invalid_ = false;

    if (artificial_ecm_treatment_)
    {
      ///////////////////////////////////////////////////////////////////////////////////
      // remove the effect of the ecm mechanics in the artificial domain
      // 1) we do not assemble to the pseudo boundary from the artificial domain.
      // 2) we fix the the artificial ECM nodes.
      // 3) we solve the real ECM part.
      // 4) we assemble as ususal.
      // 5) we fix the real ECM nodes.
      // 6) we perform a relaxation solve to solve for the artifical displacements.
      ///////////////////////////////////////////////////////////////////////////////////
      Teuchos::RCP<Epetra_Map> additional_dbc_dofmap_ecm = Teuchos::null;

      // apply displacement dbc to artifical ECM nodes
      ApplyDirichletToArtificialECM(poroscatra_subproblem_->PoroField()->StructureField(),
          backgroundstructuredis_, "PoroCoupling", additional_dbc_dofmap_ecm, 3,
          poroscatra_subproblem_->PoroField()->StructureField()->Dispnp());

      // rebuild the monolithic poro dbc map
      poroscatra_subproblem_->PoroField()->BuildCombinedDBCMap();

      // we do not assemble to pseudo-boundary during this solve
      exchange_manager_->SetPseudoBoundarySwitch(true);

      // solve poro
      poroscatra_subproblem_->Solve();

      // remove dbc from artifical ECM
      RemoveDirichlet(
          additional_dbc_dofmap_ecm, poroscatra_subproblem_->PoroField()->StructureField());

      // we do assemble to pseudo-boundary during following relaxation solve
      exchange_manager_->SetPseudoBoundarySwitch(false);

      // reuse map
      additional_dbc_dofmap_ecm = Teuchos::null;

      // apply displacement dbc to real ECM nodes
      ApplyDirichletToRealECM(poroscatra_subproblem_->PoroField()->StructureField(),
          backgroundstructuredis_, "PoroCoupling", additional_dbc_dofmap_ecm, 3,
          poroscatra_subproblem_->PoroField()->StructureField()->Dispnp());

      // rebuild the monolithic poro dbc map
      poroscatra_subproblem_->PoroField()->BuildCombinedDBCMap();

      // relaxation solve of artificial domain
      poroscatra_subproblem_->Solve();

      // remove dbc from artifical ECM
      RemoveDirichlet(
          additional_dbc_dofmap_ecm, poroscatra_subproblem_->PoroField()->StructureField());
    }
    else
      poroscatra_subproblem_->Solve();

    // in general after solving the ECM, we need to update the
    // immersed information. In case of Dirichlet coupling we
    // assume the adhesion nodes nodes not to change.
    if (!dirichletcoupling_) immersed_information_invalid_ = true;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedAdhesionTraction::ImmersedOp(
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
    // declare map for immersed adhesion dofs
    Teuchos::RCP<Epetra_Map> dirichmap_immersed;

    if (dirichletcoupling_)
    {
      // prescribe dirichlet values at cell adhesion nodes
      CalcAdhesionDisplacements();

      Teuchos::RCP<const Epetra_Map> dirichmap_original =
          cellstructure_->GetDBCMapExtractor()->CondMap();

      // build map of immersed adhesion nodes to be subjected to dirichlet conditions
      BuildImmersedDirichMap(
          immerseddis_, dirichmap_original, dirichmap_immersed, adh_nod_backgrd_ele_mapping_);
      // add adhesion dofs to dbc map
      cellstructure_->AddDirichDofs(dirichmap_immersed);
      // write dirichlet values to displacement vector
      DoImmersedDirichletCond(cellstructure_->WriteAccessDispnp(), cell_adhesion_disp_,
          cellstructure_->GetDBCMapExtractor()->CondMap());
    }

    // set state in nox and gobal state object
    cellstructure_->SetState(cellstructure_->WriteAccessDispnp());

    // solve cell
    if (Comm().MyPID() == 0)
    {
      std::cout << "\n****************************************\n          ADHESION "
                   "FORMATION\n****************************************\n";
    }
    cellscatra_subproblem_->OuterLoop();

    // remove dirichlet conditions from current adhesion nodes (set may change)
    if (dirichletcoupling_) cellstructure_->RemoveDirichDofs(dirichmap_immersed);

    return cellstructure_->ExtractImmersedInterfaceDispnp();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::PrepareBackgroundOp()
{
  immerseddis_->SetState(0, "displacement", cellstructure_->Dispnp());
  backgroundfluiddis_->SetState(0, "dispnp", poroscatra_subproblem_->FluidField()->Dispnp());

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
    immerseddis_->Dof(0, node, 0, dofstoextract);
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

  if (myrank_ == 0) std::cout << "\nPrepareBackgroundOp returns: " << std::endl;
  Comm().Barrier();
  if (curr_subset_of_backgrounddis_.empty() == false)
    std::cout << "  " << curr_subset_of_backgrounddis_.begin()->second.size()
              << " background elements on Proc " << Comm().MyPID() << std::endl;
  Comm().Barrier();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  immerseddis_->ClearState();
  immersedscatradis_->ClearState();
  backgroundfluiddis_->ClearState();
  backgroundstructuredis_->ClearState();

  // reset adhesion fixpoints (only needed for penalty method)
  if (!dirichletcoupling_ and not determine_adhesion_information_) UpdateAdhesionInformation();

  // set state to be available in bond reactions
  if (!dirichletcoupling_)
    immersedscatradis_->SetState(1, "AdhesionFixpoints", cell_adhesion_nod_coords_);

  if (myrank_ == 0) std::cout << "Cell-SSI Predictor: " << std::endl;
  cellscatra_subproblem_->PrepareTimeStep(false);
  if (myrank_ == 0) std::cout << "Poro-Scatra Predictor: " << std::endl;
  poroscatra_subproblem_->PrepareTimeStep(false);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedAdhesionTraction::InitialGuess()
{
  if (dirichletcoupling_)
    return cellstructure_->Freact();
  else
    return cell_adhesion_forces_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::UpdateAdhesionInformation()
{
  if (myrank_ == 0)
  {
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
    std::cout << "###   Update Adhesion Fixpoints ...                  " << std::endl;
  }

  const int num_boundspecies =
      globalproblem_->CellMigrationParams().sublist("ADHESION MODULE").get<int>("NUM_BOUNDSPECIES");

  // tolerance for geometric operations
  const double tol = 1e-11;

  // declare element pointers
  DRT::Element* ele = NULL;  //!< pointer to background structure ecm element
  DRT::Element* iele =
      NULL;  //!< pointer to background fluid element (carries immersed information)

  // declarte node pointer
  DRT::Node* adhesion_node_ptr = NULL;  //!< auxiliary pointer to adhesion node

  // necessary variables
  LINALG::Matrix<3, 1> xi;         //!< parameter space coordinate of anode in background element
  LINALG::Matrix<8, 1> shapefcts;  //!< shapefunctions must be evaluated at xi for force spreading
  std::vector<double> anode_coord(3, 0.0);  //!< current coordinates of anode
  std::vector<double>
      myvalues;  //!< aux. variable for extraction of element values from state vectors
  std::vector<double>
      myvalues_conc;  //!< aux. variable for extraction of element values from state vectors

  // local declarations and definitions for round robin loop
  const int numproc = numproc_;                           //!< total number of processors
  const int myrank = myrank_;                             //!< my processor id
  const int torank = (myrank + 1) % numproc;              //!< sends to
  const int fromrank = (myrank + numproc - 1) % numproc;  //!< recieves from
  int origin = myrank;                                    //!< owner of given point
  DRT::Exporter exporter(Comm());                         //!< construct eporter object

  // auxiliary counter for loop execution
  int anode_count = -1;

  // set states on participating discretizations
  immerseddis_->SetState(0, "displacement", cellstructure_->Dispnp());
  immerseddis_->SetState(1, "Phinp", cellscatra_subproblem_->ScaTraField()->ScaTraField()->Phinp());
  backgroundstructuredis_->SetState(0, "displacement", porostructure_->Dispnp());

  // get condition which marks adhesion nodes
  DRT::Condition* condition = immerseddis_->GetCondition("CellFocalAdhesion");
  const std::vector<int>* adhesion_node_ptrs = condition->Nodes();
  int my_adhesion_node_ptrsize = adhesion_node_ptrs->size();

  // get max size over all procs
  int max_adhesion_node_ptr_size_over_procs = -1;
  Comm().MaxAll(&my_adhesion_node_ptrsize, &max_adhesion_node_ptr_size_over_procs, 1);

  Teuchos::RCP<const Epetra_Vector> ecmstate = backgroundstructuredis_->GetState("displacement");
  if (ecmstate == Teuchos::null)
    dserror("Could not get state displacement from background structure");
  Teuchos::RCP<const Epetra_Vector> cellstate = immerseddis_->GetState("displacement");
  if (cellstate == Teuchos::null) dserror("Could not get state displacement from cell structure");
  Teuchos::RCP<const Epetra_Vector> cell_conc = immerseddis_->GetState(1, "Phinp");
  if (cell_conc == Teuchos::null) dserror("Could not get state Phinp from cell structure");


  // loop over all adhesion nodes
  for (int anode = 0; anode < max_adhesion_node_ptr_size_over_procs; anode++)
  {
    // has this adhesion node been matched on any proc?
    int matched = 0;
    // result matching id of b
    int matching_background_ele_id = -1;
    // determine new anchor position for this adhesion node
    int determine_new_position = 0;

    // if no more own adhesion nodes are left to match on calling proc,
    // send last node around until end. This way calling proc enters
    // round robin loop below and communication can take place.
    if (anode > (my_adhesion_node_ptrsize - 1))
    {
      anode = my_adhesion_node_ptrsize - 1;
      determine_new_position = 0;
      matched = 1;
    }
    else
      anode_count = anode;

    // get adhesion node id
    int anodeid = adhesion_node_ptrs->at(anode_count);
    // get node pointer from immersed discretization
    adhesion_node_ptr = immerseddis_->gNode(anodeid);

    // set parameter space coordinates to arbitrary not matching values
    xi(0) = 2.0;
    xi(1) = 2.0;
    xi(2) = 2.0;
    // get coordinates of adhesion node
    const double* X = adhesion_node_ptr->X();

    // get displacements of cell
    DRT::Element** adjacent_elements = adhesion_node_ptr->Elements();
    DRT::Element::LocationArray la(2);
    adjacent_elements[0]->LocationVector(*immerseddis_, la, false);
    // extract local values of the global vectors
    myvalues.resize(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*cellstate, myvalues, la[0].lm_);

    // get concentrations
    // extract local values of the global vectors
    myvalues_conc.resize(la[1].lm_.size());
    DRT::UTILS::ExtractMyValues(*cell_conc, myvalues_conc, la[1].lm_);
    double adhesioneleconc[la[1].lm_.size() * 4];
    const int numdofpernode = la[1].lm_.size() / 4;
    for (int node = 0; node < numdofpernode; ++node)
      for (int dof = 0; dof < 2; ++dof)
        adhesioneleconc[node * 2 + dof] = myvalues_conc[node * 2 + dof];

    // determine which node of element is anode
    int locid = -1;
    DRT::Node** nodes = adjacent_elements[0]->Nodes();
    for (int i = 0; i < adjacent_elements[0]->NumNode(); ++i)
    {
      if (nodes[i]->Id() == anodeid)
      {
        locid = i;
        break;
      }
    }
    if (locid == -1) dserror("could not get local index of adhesion node in element");

    // concentration of bound surface proteins
    double anode_conc = 0.0;
    for (int ii = 0; ii < num_boundspecies; ++ii)
      anode_conc += adhesioneleconc[locid * numdofpernode + ii];
    // if anode_conc > reset_conc_ --> adhesion site stays at previous site
    // otherwise we have to reset the adhesion position in the background medium
    if (anode_conc > reset_conc_ and adh_nod_param_coords_in_backgrd_ele_.find(anodeid) !=
                                         adh_nod_param_coords_in_backgrd_ele_.end())
      determine_new_position = 0;
    else if (anode_conc <= reset_conc_ and adh_nod_param_coords_in_backgrd_ele_.find(anodeid) !=
                                               adh_nod_param_coords_in_backgrd_ele_.end())
      determine_new_position = 1;
    else if (anode_conc > reset_conc_ and
             adh_nod_param_coords_in_backgrd_ele_.find(anodeid) ==
                 adh_nod_param_coords_in_backgrd_ele_.end() and
             adhesion_node_ptr->Owner() == myrank_)
      determine_new_position = 1;

    // fill vector with current coordinate
    anode_coord[0] = X[0] + myvalues[locid * 3 + 0];
    anode_coord[1] = X[1] + myvalues[locid * 3 + 1];
    anode_coord[2] = X[2] + myvalues[locid * 3 + 2];


    ///////////////////////////////////////////////////////////////////////////////////////////
    ////
    ////     round robin loop
    ////
    //// find ecm element in which adhesion node is immersed. this ecm element is not guaranteed
    //// to be on the same proc as the adhesion node.
    ////
    //// in this round-robin loop, every proc sends the adhesion node coordinates to the next
    //// proc. Each proc tries to match the recieved coordinates to a background element.
    ////
    //// when the loop ends, all information are back on the proc who requested the data
    //// originally and for every point the data should be stored in rdata.
    ////
    //// --> loop from 0 to numproc instead of (numproc-1) is chosen intentionally.
    ////
    ////////////////////////////////////////////////////////////////////////////////////////////
    // round robin loop
    for (int irobin = 0; irobin < numproc; ++irobin)
    {
      std::vector<char> sdata;
      std::vector<char> rdata;

      // each proc searches in its own subsetof backgrd. elements for a match with the adhesion node
      for (std::map<int, std::set<int>>::const_iterator closele =
               curr_subset_of_backgrounddis_.begin();
           closele != curr_subset_of_backgrounddis_.end(); closele++)
      {
        if (not determine_new_position or matched) break;

        // search in each background element
        for (std::set<int>::const_iterator eleIter = (closele->second).begin();
             eleIter != (closele->second).end(); eleIter++)
        {
          if (not determine_new_position or matched) break;

          ele = backgroundstructuredis_->gElement(*eleIter);  //< background solid phase element
          iele = backgroundfluiddis_->gElement(
              *eleIter);  //< background fluid element with same id has all immersed information

          DRT::ELEMENTS::FluidImmersedBase* immersedele =
              dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(iele);
          if (immersedele == NULL)
            dserror("dynamic cast from DRT::Element* to DRT::ELEMENTS::FluidImmersedBase* failed");

          // if we need to determine the local coordinates of the adhesion nodes
          // we only try to match owned adhesion nodes (row layout of mappings)
          if (determine_new_position)
          {
            bool converged = false;
            double residual = -1234.0;
            matching_background_ele_id = -1234;

            // get background displacements
            DRT::Element::LocationArray la(1);
            ele->LocationVector(*backgroundstructuredis_, la, false);
            // extract local values of the global vectors
            myvalues.resize(la[0].lm_.size());
            DRT::UTILS::ExtractMyValues(*ecmstate, myvalues, la[0].lm_);
            double sourceeledisp[24];
            for (int node = 0; node < 8; ++node)
              for (int dof = 0; dof < 3; ++dof)
                sourceeledisp[node * 3 + dof] = myvalues[node * 4 + dof];

            // node 1  and node 7 coords of current backgrd. element (diagonal points)
            const double* X1 = ele->Nodes()[1]->X();
            double x1[3];
            x1[0] = X1[0] + sourceeledisp[1 * 3 + 0];
            x1[1] = X1[1] + sourceeledisp[1 * 3 + 1];
            x1[2] = X1[2] + sourceeledisp[1 * 3 + 2];
            const double* X7 = ele->Nodes()[7]->X();
            double diagonal =
                sqrt(pow(X1[0] - X7[0], 2) + pow(X1[1] - X7[1], 2) + pow(X1[2] - X7[2], 2));

            // calc distance of current anode to arbitrary node (e.g. node 1) of curr source element
            double distance = sqrt(pow(x1[0] - anode_coord[0], 2) + pow(x1[1] - anode_coord[1], 2) +
                                   pow(x1[2] - anode_coord[2], 2));

            // get parameter space coords xi in source element of global point anode
            // NOTE: if the anode is very far away from the source element ele
            //       it is unnecessary to jump into this method and invoke a newton iteration.
            // Therefore: only call GlobalToCurrentLocal if distance is smaller than
            // factor*characteristic element length
            if (distance < 2.5 * diagonal)
            {
              MORTAR::UTILS::GlobalToCurrentLocal<DRT::Element::hex8>(
                  *ele, &sourceeledisp[0], &anode_coord[0], &xi(0), converged, residual);

              if (converged == false)
              {
                std::cout << "Warning! GlobalToCurrentLocal did not converge for adhesion node "
                          << anodeid << ". Res=" << residual << std::endl;
                xi(0) = 2.0;
                xi(1) = 2.0;
                xi(2) = 2.0;
              }
            }
            else
            {
              xi(0) = 2.0;
              xi(1) = 2.0;
              xi(2) = 2.0;
            }

            // anode lies in element ele so we found a match.
            // we only match adhesion node if this proc owns the background element
            if (abs(xi(0)) < (1.0 + tol) and abs(xi(1)) < (1.0 + tol) and
                abs(xi(2)) < (1.0 + tol) and ele->Owner() == myrank_)
            {
              matched = 1;
              matching_background_ele_id = ele->Id();
              break;
            }  // if match
          }    // if determine_new_position

        }  // loop over all background elements
      }    // loop over curr_subset_of_backgrounddis_


      if (numproc > 1)
      {
        // ---- pack data for sending -----
        {
          DRT::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(anodeid);
            data.AddtoPack(matching_background_ele_id);
            data.AddtoPack(origin);
            data.AddtoPack(determine_new_position);
            data.AddtoPack(matched);

            // pack vectors
            for (int dim = 0; dim < 3; ++dim)
            {
              data.AddtoPack(anode_coord[dim]);
              data.AddtoPack(xi(dim));
            }
          }

          std::swap(sdata, data());
        }
        //---------------------------------

        // ---- send ----
        MPI_Request request;
        exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from, tag, rdata, length);
        if (tag != 1234 or from != fromrank)
          dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank,
              myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;

        anodeid = DRT::ParObject::ExtractInt(position, rdata);
        matching_background_ele_id = DRT::ParObject::ExtractInt(position, rdata);
        origin = DRT::ParObject::ExtractInt(position, rdata);
        determine_new_position = DRT::ParObject::ExtractInt(position, rdata);
        matched = DRT::ParObject::ExtractInt(position, rdata);

        for (int dim = 0; dim < 3; ++dim)
        {
          anode_coord[dim] = DRT::ParObject::ExtractDouble(position, rdata);
          xi(dim) = DRT::ParObject::ExtractDouble(position, rdata);
        }

        // wait for all communication to finish
        exporter.Wait(request);
        Comm().Barrier();
      }  // communicate data only if num procs > 1

    }  // round robin loop


    // now all data should be back at originating processor
    if (origin != myrank_)
      dserror("After round-robin loop the delivered data should arrive at the origin again.");
    if (determine_new_position and (not matched))
      dserror("Adhesion node with id %i could not be matched on any proc. %f<%f , %i", anodeid,
          anode_conc, reset_conc_, determine_new_position);

    // adhesion node was matched to a background element and xi coordinate was communicated.
    // 1) we replace the adhesion information for this node.
    // 2) we write the current global coordinates of the adhesion node to cell_adhesion_nod_coords_.
    if (matched and determine_new_position and adhesion_node_ptr->Owner() == myrank_)
    {
      // delete old pair cell adhesion node id -> xi in backgroundele in map
      adh_nod_param_coords_in_backgrd_ele_.erase(anodeid);
      // delete old pair adhesion node id -> backgrdele id in map
      adh_nod_backgrd_ele_mapping_.erase(anodeid);
      // write new pair cell adhesion node id -> xi in backgroundele in map
      adh_nod_param_coords_in_backgrd_ele_.insert(
          std::pair<int, LINALG::Matrix<3, 1>>(anodeid, xi));
      // write new pair adhesion node id -> backgrdele id in map
      adh_nod_backgrd_ele_mapping_.insert(std::pair<int, int>(anodeid, ele->Id()));

      // write adhesion fixpoint coordinates in cell_adhesion_nod_coords_

      // get dofs of cell adhesion node
      std::vector<int> dofs = immerseddis_->Dof(0, immerseddis_->gNode(anodeid));
      if (dofs.size() != 3) dserror("dofs=3 expected. dofs=%d instead", dofs.size());

      // write displacmement into global vector
      for (int dof = 0; dof < 3; dof++)
      {
        // get lid of current gid in 'dofs'
        int lid = immerseddis_->DofRowMap()->LID(dofs[dof]);
        int err = 0;
        if (lid > -1) err = cell_adhesion_nod_coords_->ReplaceMyValue(lid, 0, anode_coord[dof]);
        if (err != 0) dserror("ReplaceMyValue returned err=%d", err);
      }  // loop over all dofs of adhesion anode

    }  // adhesion node is matched

  }  // loop over all nodes with adhesion condition 'anode'


  // set state to be available in bond reactions
  immersedscatradis_->SetState(1, "AdhesionFixpoints", cell_adhesion_nod_coords_);

  Comm().Barrier();
  std::cout << "###   PROC " << myrank_ << ": Size of 'adh_nod_param_coords_in_backgrd_ele_' = "
            << adh_nod_param_coords_in_backgrd_ele_.size() << std::endl;
  Comm().Barrier();
  std::cout << "###   PROC " << myrank_ << ": Size of 'adh_nod_backgrd_ele_mapping_'         = "
            << adh_nod_backgrd_ele_mapping_.size() << std::endl;
  Comm().Barrier();

  double norm = -1234.0;
  cell_adhesion_nod_coords_->Norm1(&norm);

  if (myrank_ == 0)
  {
    std::cout << "###   L1-Norm of adhesion node coords: " << std::setprecision(11) << norm
              << std::endl;
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
  }

  return;
}  // UpdateAdhesionInformation()


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::DistributeAdhesionForce(
    const Teuchos::RCP<const Epetra_Vector>& force_to_distribute)
{
  if (myrank_ == 0)
  {
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
    std::cout << "###   Spread adhesion forces onto ecm ...                  " << std::endl;
  }

  // safety check
  if (immersed_information_invalid_) dserror("immersed information should be updated first");

  // reinitialize map to insert new adhesion node -> background element id pairs
  if (determine_adhesion_information_)
  {
    adh_nod_backgrd_ele_mapping_.clear();
    adh_nod_param_coords_in_backgrd_ele_.clear();
  }

  // tolerance for geometric operations
  const double tol = 1e-11;

  // initialize vector norms
  double adhesionforcenorm = -1234.0;
  double freactnorm = -1234.0;

  // declare element pointers
  DRT::Element* ele = NULL;  //!< pointer to background structure ecm element
  DRT::Element* iele =
      NULL;  //!< pointer to background fluid element (carries immersed information)

  // declarte node pointer
  DRT::Node* adhesion_node = NULL;  //!< auxiliary pointer to adhesion node

  // necessary variables
  LINALG::Matrix<3, 1> xi;         //!< parameter space coordinate of anode in background element
  LINALG::Matrix<8, 1> shapefcts;  //!< shapefunctions must be evaluated at xi for force spreading
  std::vector<double> anode_freact(3, 0.0);  //!< reaction force at adhesion node
  std::vector<double> anode_coord(3, 0.0);   //!< current coordinates of anode
  std::vector<double>
      myvalues;  //!< aux. variable for extraction of element values from state vectors

  // local declarations and definitions for round robin loop
  const int numproc = numproc_;                           //!< total number of processors
  const int myrank = myrank_;                             //!< my processor id
  const int torank = (myrank + 1) % numproc;              //!< sends to
  const int fromrank = (myrank + numproc - 1) % numproc;  //!< recieves from
  int origin = myrank;                                    //!< owner of given point
  DRT::Exporter exporter(Comm());                         //!< construct exporter object

  // auxiliary counter for loop execution
  int anode_count = -1;
  // result matching id
  int matching_id = -1;
  // current adhesion node already matched?
  int finished = 0;

  // aux variables to be communicated to add contribution of
  // column background nodes on each proc:
  // an adhesion node my be owned on proc X and the covering
  // backgrd. ele on proc Y. some of the nodes of this ele may
  // however be owned by proc Z. Since the vector ecm_adhesion_force
  // has row layout, we can add entries only on procs owning the node.
  // So we have to communicate the force from proc Y to proc Z and add
  // the force contribution on proc Z.
  std::vector<int> column_node_id;
  std::vector<LINALG::Matrix<3, 1>> column_node_force;

  // set states on participating discretizations
  immerseddis_->SetState(0, "displacement", cellstructure_->Dispnp());
  backgroundstructuredis_->SetState(0, "displacement", porostructure_->Dispnp());

  // get reaction forces at cell surface
  double* freact_values = force_to_distribute->Values();
  force_to_distribute->Norm1(&freactnorm);
  if (myrank_ == 0)
    std::cout << "###   L1-Norm of Cell Reaction Force Vector: " << std::setprecision(11)
              << freactnorm << std::endl;

  // get condition which marks adhesion nodes
  DRT::Condition* condition = immerseddis_->GetCondition("CellFocalAdhesion");
  const std::vector<int>* adhesion_nodes = condition->Nodes();
  int my_adhesion_nodesize = adhesion_nodes->size();

  // get max size ove all procs
  int max_adhesion_node_size_over_procs = -1;
  Comm().MaxAll(&my_adhesion_nodesize, &max_adhesion_node_size_over_procs, 1);

  Teuchos::RCP<const Epetra_Vector> ecmstate = backgroundstructuredis_->GetState("displacement");
  if (ecmstate == Teuchos::null)
    dserror("Could not get state displacement from background structure");

  Teuchos::RCP<const Epetra_Vector> cellstate = immerseddis_->GetState("displacement");
  if (cellstate == Teuchos::null) dserror("Could not get state displacement from cell structure");

  // if immersed information are invalid, we have to match adhesion nodes and
  // background elements again.
  // loop over all adhesion nodes
  // each proc needs to enter here until the proc with the most
  // adhesion nodes has finished, because each proc needs to run
  // into the round robin loop to make it work.
  // todo cell is ghosted redundantly. max_adhesion_node_size_over_procs
  // equals my_adhesion_nodesize anyway. We should only loop row nodes!
  // this requires most likely further changes to run.
  for (int anode = 0; anode < max_adhesion_node_size_over_procs; anode++)
  {
    // reset matching id
    matching_id = -1;
    // reset matched flag
    finished = 0;

    // if no more own adhesion nodes are left to match on calling proc,
    // send last node around until end. This way calling proc enters
    // round robin loop below and communication can take place.
    if (anode > (my_adhesion_nodesize - 1))
    {
      anode = my_adhesion_nodesize - 1;
      finished = 1;
    }
    else
      anode_count = anode;

    // get adhesion node id
    int anodeid = adhesion_nodes->at(anode_count);

    // get node pointer from immersed discretization
    adhesion_node = immerseddis_->gNode(anodeid);


    if (not determine_adhesion_information_ and adhesion_node->Owner() == myrank_)
    {
      xi = adh_nod_param_coords_in_backgrd_ele_.at(anodeid);
      matching_id = adh_nod_backgrd_ele_mapping_.at(anodeid);
      if (matching_id < 0) dserror("invalid id for matching background element in map");
    }
    else
    {
      // set parameter space coordinates to arbitrary not matching values
      xi(0) = 2.0;
      xi(1) = 2.0;
      xi(2) = 2.0;
    }

    // get coordinates of adhesion node
    const double* X = adhesion_node->X();

    // get nodal adhesion force
    std::vector<int> dofs = immerseddis_->Dof(0, adhesion_node);
    if (dofs.size() != 3) dserror("dofs=3 expected. dofs=%d instead", dofs.size());

    for (int dim = 0; dim < (3 - isPseudo2D_); dim++)
    {
      int doflid = force_to_distribute->Map().LID(dofs[dim]);
      anode_freact[dim] = freact_values[doflid];
    }  // dof loop

    // get displacements of cell
    DRT::Element** adjacent_elements = adhesion_node->Elements();
    DRT::Element::LocationArray la(1);
    adjacent_elements[0]->LocationVector(*immerseddis_, la, false);
    // extract local values of the global vectors
    myvalues.resize(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*cellstate, myvalues, la[0].lm_);

    // determine which node of element is anode
    int locid = -1;
    DRT::Node** nodes = adjacent_elements[0]->Nodes();
    for (int i = 0; i < adjacent_elements[0]->NumNode(); ++i)
    {
      if (nodes[i]->Id() == anodeid)
      {
        locid = i;
        break;
      }
    }

    if (locid == -1) dserror("could not get local index of adhesion node in element");

    // fill vector with current coordinate
    anode_coord[0] = X[0] + myvalues[locid * 3 + 0];
    anode_coord[1] = X[1] + myvalues[locid * 3 + 1];
    anode_coord[2] = X[2] + myvalues[locid * 3 + 2];

    ///////////////////////////////////////////////////////////////////////////////////////////
    ////
    ////     round robin loop
    ////
    //// find ecm element in which adhesion node is immersed. this ecm element is not guaranteed
    //// to be on the same proc as the adhesion node.
    ////
    //// in this round-robin loop, every proc sends the adhesion node coordinates to the next
    //// proc. Each proc tries to match the recieved coordinates to a background element.
    ////
    //// when the loop ends, all information are back on the proc who requested the data
    //// originally and for every point the data should be stored in rdata.
    ////
    //// --> loop from 0 to numproc instead of (numproc-1) is chosen intentionally.
    ////
    ////////////////////////////////////////////////////////////////////////////////////////////
    // round robin loop
    for (int irobin = 0; irobin < numproc; ++irobin)
    {
      std::vector<char> sdata;
      std::vector<char> rdata;

      bool mymatch = false;

      // each proc searches in its own subsetof backgrd. elements for a match with the adhesion node
      for (std::map<int, std::set<int>>::const_iterator closele =
               curr_subset_of_backgrounddis_.begin();
           closele != curr_subset_of_backgrounddis_.end(); closele++)
      {
        // if the adhesion node is already matched, we can skip the following
        if (finished > 0) break;

        // search in each background element
        for (std::set<int>::const_iterator eleIter = (closele->second).begin();
             eleIter != (closele->second).end(); eleIter++)
        {
          // if the adhesion node is already matched, we can skip the following
          if (finished > 0) break;

          ele = backgroundstructuredis_->gElement(*eleIter);  //< background solid phase element
          iele = backgroundfluiddis_->gElement(
              *eleIter);  //< background fluid element with same id has all immersed information

          DRT::ELEMENTS::FluidImmersedBase* immersedele =
              dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(iele);
          if (immersedele == NULL)
            dserror("dynamic cast from DRT::Element* to DRT::ELEMENTS::FluidImmersedBase* failed");

          // if we need to determine the local coordinates of the adhesion nodes
          // match is only possible if background element is cut by immersed boundary
          // we only try to matched owned adhesion nodes (row layout of mappings)
          if (determine_adhesion_information_ and immersedele->IsBoundaryImmersed())
          {
            bool converged = false;
            double residual = -1234.0;

            // reset matching id
            matching_id = -1234;

            // get background displacements
            DRT::Element::LocationArray la(1);
            ele->LocationVector(*backgroundstructuredis_, la, false);
            // extract local values of the global vectors
            myvalues.resize(la[0].lm_.size());
            DRT::UTILS::ExtractMyValues(*ecmstate, myvalues, la[0].lm_);
            double sourceeledisp[24];
            for (int node = 0; node < 8; ++node)
              for (int dof = 0; dof < 3; ++dof)
                sourceeledisp[node * 3 + dof] = myvalues[node * 4 + dof];

            // node 1  and node 7 coords of current backgrd. element (diagonal points)
            const double* X1 = ele->Nodes()[1]->X();
            double x1[3];
            x1[0] = X1[0] + sourceeledisp[1 * 3 + 0];
            x1[1] = X1[1] + sourceeledisp[1 * 3 + 1];
            x1[2] = X1[2] + sourceeledisp[1 * 3 + 2];
            const double* X7 = ele->Nodes()[7]->X();
            double diagonal =
                sqrt(pow(X1[0] - X7[0], 2) + pow(X1[1] - X7[1], 2) + pow(X1[2] - X7[2], 2));

            // calc distance of current anode to arbitrary node (e.g. node 1) of curr source element
            double distance = sqrt(pow(x1[0] - anode_coord[0], 2) + pow(x1[1] - anode_coord[1], 2) +
                                   pow(x1[2] - anode_coord[2], 2));

            // get parameter space coords xi in source element of global point anode
            // NOTE: if the anode is very far away from the source element ele
            //       it is unnecessary to jump into this method and invoke a newton iteration.
            // Therefore: only call GlobalToCurrentLocal if distance is smaller than
            // factor*characteristic element length
            if (distance < 2.5 * diagonal)
            {
              MORTAR::UTILS::GlobalToCurrentLocal<DRT::Element::hex8>(
                  *ele, &sourceeledisp[0], &anode_coord[0], &xi(0), converged, residual);

              if (converged == false)
              {
                std::cout << "Warning! GlobalToCurrentLocal did not converge for adhesion node "
                          << anodeid << ". Res=" << residual << std::endl;
                xi(0) = 2.0;
                xi(1) = 2.0;
                xi(2) = 2.0;
              }
            }
            else
            {
              xi(0) = 2.0;
              xi(1) = 2.0;
              xi(2) = 2.0;
            }

            // anode lies in element ele so we found a match.
            // we only match adhesion node if this proc owns the background element
            if (abs(xi(0)) < (1.0 + tol) and abs(xi(1)) < (1.0 + tol) and
                abs(xi(2)) < (1.0 + tol) and ele->Owner() == myrank_)
            {
              mymatch = true;
              matching_id = ele->Id();
            }
          }
          else
          {
            // xi and matching_id already set above
            // we just check if ele has matching_id
            if (ele->Id() == matching_id and ele->Owner() == myrank_)
            {
              mymatch = true;
            }
          }

          // we spread the adhesion force of the adhesion node to
          // the nodes of the matching background element.
          if (mymatch)
          {
            // spread force to nodes of ecm ele
            DRT::Node** spreadnodes = ele->Nodes();
            // evaluate shapefcts of ele at point xi
            DRT::UTILS::shape_function<DRT::Element::hex8>(xi, shapefcts);

            for (int snode = 0; snode < (int)ele->NumNode(); snode++)
            {
              bool owned = false;
              if (spreadnodes[snode]->Owner() == myrank_) owned = true;

              std::vector<int> sdofs = backgroundstructuredis_->Dof(0, spreadnodes[snode]);
              if (sdofs.size() != 4)
                dserror("dofs=4 expected. dofs=%d instead",
                    sdofs.size());  // 4 dofs per node in porostructure

              if (owned)
              {
                // write entry into ecm_adhesion_forces_
                for (int dim = 0; dim < (3 - isPseudo2D_); dim++)
                {
                  int error = ecm_adhesion_forces_->SumIntoGlobalValue(
                      sdofs[dim], 0, 0, shapefcts(snode) * anode_freact[dim]);
                  if (error != 0) dserror("SumIntoGlobalValue returned err=%d", error);

                }  // dim loop
              }    // if node owned
              else
              {
                LINALG::Matrix<3, 1> tempforce;
                tempforce(0) = shapefcts(snode) * anode_freact[0];
                tempforce(1) = shapefcts(snode) * anode_freact[1];
                tempforce(2) = shapefcts(snode) * anode_freact[2];

                column_node_id.push_back(spreadnodes[snode]->Id());
                column_node_force.push_back(tempforce);
              }
            }  // node loop

            finished = true;
          }  // found match

        }  // loop over all background elements
      }    // loop over curr_subset_of_backgrounddis_

      if (numproc > 1)
      {
        // ---- pack data for sending -----
        {
          DRT::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(anodeid);
            data.AddtoPack(finished);
            data.AddtoPack(matching_id);
            data.AddtoPack(origin);

            // pack vectors
            for (int dim = 0; dim < 3; ++dim)
            {
              data.AddtoPack(anode_freact[dim]);
              data.AddtoPack(anode_coord[dim]);
              data.AddtoPack(xi(dim));
            }
          }

          std::swap(sdata, data());
        }
        //---------------------------------

        // ---- send ----
        MPI_Request request;
        exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from, tag, rdata, length);
        if (tag != 1234 or from != fromrank)
          dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank,
              myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;

        anodeid = DRT::ParObject::ExtractInt(position, rdata);
        finished = DRT::ParObject::ExtractInt(position, rdata);
        matching_id = DRT::ParObject::ExtractInt(position, rdata);
        origin = DRT::ParObject::ExtractInt(position, rdata);

        for (int dim = 0; dim < 3; ++dim)
        {
          anode_freact[dim] = DRT::ParObject::ExtractDouble(position, rdata);
          anode_coord[dim] = DRT::ParObject::ExtractDouble(position, rdata);
          xi(dim) = DRT::ParObject::ExtractDouble(position, rdata);
        }

        // wait for all communication to finish
        exporter.Wait(request);
        Comm().Barrier();
      }  // communicate data only if num procs > 1
    }    // round robin loop

    // now all data should be back at originating processor
    if (origin != myrank_)
      dserror("After round-robin loop the delivered data should arrive at the origin again.");
    // we only have to match row nodes
    if (not finished and adhesion_node->Owner() == myrank_)
      dserror(
          "Adhesion node %i could not be matched to a background element on any proc.\n"
          "This node has initial coordinates [%f %f %f]",
          anodeid, adhesion_node->X()[0], adhesion_node->X()[1], adhesion_node->X()[2]);

    if (finished and adhesion_node->Owner() == myrank_ and determine_adhesion_information_)
    {
      // write pair cell adhesion node id -> xi in backgroundele in map
      adh_nod_param_coords_in_backgrd_ele_.insert(
          std::pair<int, LINALG::Matrix<3, 1>>(anodeid, xi));
      // write pair adhesion node id -> backgrdele id in map
      adh_nod_backgrd_ele_mapping_.insert(std::pair<int, int>(anodeid, matching_id));

      if (!dirichletcoupling_)
      {
        // write cordinate of adhesion node
        DRT::Node* temp_node_ptr = immerseddis_->gNode(anodeid);
        std::vector<int> dofs = immerseddis_->Dof(0, temp_node_ptr);
        for (int dof = 0; dof < 3; dof++)
        {
          // get lid of current gid in 'dofs'
          int lid = immerseddis_->DofRowMap()->LID(dofs[dof]);
          int err = 0;
          if (lid > -1) err = cell_adhesion_nod_coords_->ReplaceMyValue(lid, 0, anode_coord[dof]);
          if (err != 0) dserror("ReplaceMyValue returned err=%d", err);
        }  // loop over all dofs of adhesion anode
      }
    }

  }  // loop over all nodes with adhesion condition 'anode'


  // second round robin loop to add contributions to ecm_adhesion_forces_ to the nodes
  // associated with matched elements who did not own these nodes
  int my_size = column_node_id.size();
  int max_size_over_all_procs = -1;
  Comm().MaxAll(&my_size, &max_size_over_all_procs, 1);

  for (int i = 0; i < max_size_over_all_procs; ++i)
  {
    int node_id = -1;
    LINALG::Matrix<3, 1> node_force(true);

    if (i < my_size)
    {
      node_id = column_node_id[i];
      node_force = column_node_force[i];
    }

    // round-robin loop
    for (int irobin = 0; irobin < numproc; ++irobin)
    {
      std::vector<char> sdata;
      std::vector<char> rdata;

      // do the work ...

      // check whether element is owned by this proc
      if (backgroundstructuredis_->NodeRowMap()->LID(node_id) != -1)
      {
        DRT::Node* mynode = backgroundstructuredis_->gNode(node_id);
        std::vector<int> sdofs = backgroundstructuredis_->Dof(0, mynode);
        if (sdofs.size() != 4)
          dserror("dofs=4 expected. dofs=%d instead",
              sdofs.size());  // 4 dofs per node in porostructure

        // write entry into ecm_adhesion_forces_
        for (int dim = 0; dim < (3 - isPseudo2D_); dim++)
        {
          int error = ecm_adhesion_forces_->SumIntoGlobalValue(sdofs[dim], 0, 0, node_force(dim));
          if (error != 0) dserror("SumIntoGlobalValue returned err=%d", error);
        }  // dim loop
      }    // if proc owns node


      // do the communication ...
      if (numproc > 1)
      {
        // ---- pack data for sending -----
        {
          DRT::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(origin);
            data.AddtoPack(node_id);

            for (int dim = 0; dim < 3; ++dim) data.AddtoPack(node_force(dim));
          }

          std::swap(sdata, data());
        }
        //---------------------------------

        // ---- send ----
        MPI_Request request;
        exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from, tag, rdata, length);
        if (tag != 1234 or from != fromrank)
          dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank,
              myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;

        origin = DRT::ParObject::ExtractInt(position, rdata);
        node_id = DRT::ParObject::ExtractInt(position, rdata);

        for (int dim = 0; dim < 3; ++dim)
          node_force(dim) = DRT::ParObject::ExtractDouble(position, rdata);

        // wait for all communication to finish
        exporter.Wait(request);
        Comm().Barrier();
      }  // communicate data only if num procs > 1

    }  // irobin
  }    // do for all node ids in column_node_id

  ecm_adhesion_forces_->Norm1(&adhesionforcenorm);

  if (myrank_ == 0)
    std::cout << "###   L1-Norm of  ECM Adhesion Force Vector: " << std::setprecision(11)
              << adhesionforcenorm << std::endl;

  if (determine_adhesion_information_)
  {
    Comm().Barrier();
    std::cout << "###   PROC " << myrank_ << ": Size of 'adh_nod_param_coords_in_backgrd_ele_' = "
              << adh_nod_param_coords_in_backgrd_ele_.size() << std::endl;
    Comm().Barrier();
    std::cout << "###   PROC " << myrank_ << ": Size of 'adh_nod_backgrd_ele_mapping_'         = "
              << adh_nod_backgrd_ele_mapping_.size() << std::endl;
  }

  if (!dirichletcoupling_)
  {
    double norm = -1234.0;
    cell_adhesion_nod_coords_->Norm1(&norm);
    if (myrank_ == 0)
      std::cout << "###   L1-Norm of adhesion node coords: " << std::setprecision(11) << norm
                << std::endl;
  }

  Comm().Barrier();
  if (myrank_ == 0)
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;

  // set state to be available in bond reactions
  if (determine_adhesion_information_ and !dirichletcoupling_)
    immersedscatradis_->SetState(1, "AdhesionFixpoints", cell_adhesion_nod_coords_);

  determine_adhesion_information_ = false;

  return;
}  // DistributeAdhesionForce


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::CalcAdhesionDisplacements()
{
  double adhesiondispnorm = -1234.0;

  LINALG::Matrix<3, 1> xi;            //!< parameter space coordinate of anode in background element
  LINALG::Matrix<8, 1> shapefcts;     //!< shapefunctions must be evaluated at xi
  LINALG::Matrix<3, 1> adhesiondisp;  //!< displacement of ECM element interpolated to adhesion node

  DRT::Element* backgrdele = NULL;  //!< pointer to ECM element

  // local declarations and definitions for round robin loop
  const int numproc = numproc_;                           //!< total number of processors
  const int myrank = myrank_;                             //!< my processor id
  const int torank = (myrank + 1) % numproc;              //!< sends to
  const int fromrank = (myrank + numproc - 1) % numproc;  //!< recieves from
  // owner of given point
  int origin = myrank;
  // construct eporter object
  DRT::Exporter exporter(Comm());
  // matched flag
  int matched = 0;

  // safety check
  if (adh_nod_param_coords_in_backgrd_ele_.size() != adh_nod_backgrd_ele_mapping_.size())
    dserror("Mappings must have same number of entries.");

  // get maximum number of nodes over processors
  int my_adhesion_nodesize = adh_nod_param_coords_in_backgrd_ele_.size();
  int max_adhesionnode_size_over_procs = -1;
  Comm().MaxAll(&my_adhesion_nodesize, &max_adhesionnode_size_over_procs, 1);

  if (myrank_ == 0)
  {
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
    std::cout << "###   Apply Dirichlet to Adhering Nodes                  " << std::endl;
  }
  Comm().Barrier();
  std::cout << "###   PROC " << myrank_ << ": Size of 'adh_nod_param_coords_in_backgrd_ele_' = "
            << adh_nod_param_coords_in_backgrd_ele_.size() << std::endl;
  Comm().Barrier();
  std::cout << "###   PROC " << myrank_ << ": Size of 'adh_nod_backgrd_ele_mapping_'         = "
            << adh_nod_backgrd_ele_mapping_.size() << std::endl;
  Comm().Barrier();

  // iterator for map adhesion node id <-> parameter sape coordinate in matcing background element
  std::map<int, LINALG::Matrix<3, 1>>::iterator it = adh_nod_param_coords_in_backgrd_ele_.begin();

  // loop over all cell adhesion nodes
  for (int anode = 0; anode < max_adhesionnode_size_over_procs; ++anode)
  {
    // reset matched flag
    matched = 0;

    DRT::Node* adhesion_node = NULL;

    // get adhesion node id
    int anode_id = -1;
    xi(0) = -2.0;
    xi(1) = -2.0;
    xi(2) = -2.0;
    if (it != adh_nod_param_coords_in_backgrd_ele_.end())
    {
      // get node id
      anode_id = it->first;
      // get parameter space coordinate of cell adhesion node
      xi = it->second;
      // get adhesion node
      adhesion_node = immerseddis_->gNode(anode_id);
      if (adhesion_node->Owner() != myrank_) dserror("Only row nodes should be in map.");
    }

    // get matching background element id
    int matching_id = -1;
    if (anode_id > -1) matching_id = adh_nod_backgrd_ele_mapping_.at(anode_id);

    backgroundstructuredis_->SetState(0, "displacement", porostructure_->Dispnp());
    Teuchos::RCP<const Epetra_Vector> dispnp = backgroundstructuredis_->GetState("displacement");
    if (dispnp == Teuchos::null)
      dserror("Could not get state 'displacement' from porous structure.");

    ///////////////////////////////////////////////////////////////////////////////////////////
    ////
    ////     round robin loop
    ////
    //// interpolate ecm displacement to adhesion node. since adhesion node may be on different
    //// proc than adhesion node, we need to ask the proc having the matching background element
    //// for the displacement.
    ////
    //// --> loop from 0 to numproc instead of (numproc-1) is chosen intentionally.
    ////
    ////////////////////////////////////////////////////////////////////////////////////////////
    // round robin loop
    for (int irobin = 0; irobin < numproc; ++irobin)
    {
      std::vector<char> sdata;
      std::vector<char> rdata;

      if (not matched)
      {
        // get background element if it is on this proc
        if (backgroundstructuredis_->HaveGlobalElement(matching_id))
          backgrdele = backgroundstructuredis_->gElement(matching_id);
        else
          backgrdele = NULL;

        if (backgrdele != NULL)
        {
          if (backgrdele->Owner() == myrank_)
          {
            matched = 1;

            // 1) extract background element nodal displacements
            DRT::Element::LocationArray la(1);
            backgrdele->LocationVector(*backgroundstructuredis_, la, false);
            std::vector<double> myvalues(la[0].lm_.size());
            DRT::UTILS::ExtractMyValues(*dispnp, myvalues, la[0].lm_);


            // 2) interpolate ECM displacement to parameter space coordinate occupied by the cell
            // adhesion node prior to deformation.
            //    the adhesion node is supposed to be fixed to the same material point.

            // evaluate shapefcts of ele at point xi
            DRT::UTILS::shape_function<DRT::Element::hex8>(xi, shapefcts);

            // interpolate background ele displacement to adhesion node
            adhesiondisp.Scale(0.0);
            for (int node = 0; node < backgrdele->NumNode(); node++)
            {
              adhesiondisp(0) += shapefcts(node) * myvalues[node * 4 + 0];
              adhesiondisp(1) += shapefcts(node) * myvalues[node * 4 + 1];
              adhesiondisp(2) += shapefcts(node) * myvalues[node * 4 + 2];
            }

          }  // execute if this proc owns the current background element
        }    // execute if this proc has a background element
      }      // if not yet matched

      if (numproc > 1)
      {
        // ---- pack data for sending -----
        {
          DRT::PackBuffer data;
          data.StartPacking();
          {
            data.AddtoPack(anode_id);
            data.AddtoPack(matching_id);
            data.AddtoPack(origin);
            data.AddtoPack(matched);

            // parameter space coordinate
            for (int dim = 0; dim < 3; ++dim) data.AddtoPack(xi(dim));

            if (matched)
            {
              for (int dim = 0; dim < 3; ++dim) data.AddtoPack(adhesiondisp(dim));
            }
          }
          std::swap(sdata, data());
        }
        //---------------------------------

        // ---- send ----
        MPI_Request request;
        exporter.ISend(myrank, torank, &(sdata[0]), (int)sdata.size(), 1234, request);

        // ---- receive ----
        int length = rdata.size();
        int tag = -1;
        int from = -1;
        exporter.ReceiveAny(from, tag, rdata, length);
        if (tag != 1234 or from != fromrank)
          dserror("Received data from the wrong proc soll(%i -> %i) ist(%i -> %i)", fromrank,
              myrank, from, myrank);

        // ---- unpack data -----
        std::vector<char>::size_type position = 0;
        anode_id = DRT::ParObject::ExtractInt(position, rdata);
        matching_id = DRT::ParObject::ExtractInt(position, rdata);
        origin = DRT::ParObject::ExtractInt(position, rdata);
        matched = DRT::ParObject::ExtractInt(position, rdata);

        for (int dim = 0; dim < 3; ++dim) xi(dim) = DRT::ParObject::ExtractDouble(position, rdata);

        if (matched)
        {
          for (int dim = 0; dim < 3; ++dim)
            adhesiondisp(dim) = DRT::ParObject::ExtractDouble(position, rdata);
        }

        // wait for all communication to finish
        exporter.Wait(request);
        Comm().Barrier();
      }  // communicate data only if num procs > 1

    }  // round-robin loop

    // safety checks
    if (origin != myrank_)
      dserror("After round-robin loop the delivered data should arrive at the origin again.");
    if (not matched and anode_id > -1)
      dserror(
          "Adhesion node %i could not be matched to a background element on any proc.", anode_id);

    if (adhesion_node != NULL)
    {
      if (matched)
      {
        // 3) apply previously calculated adhesion node displacement to cell

        // get dofs of cell adhesion node
        std::vector<int> dofs = immerseddis_->Dof(0, adhesion_node);
        if (dofs.size() != 3) dserror("dofs=3 expected. dofs=%d instead", dofs.size());

        bool owned = adhesion_node->Owner() == myrank_;

        if (owned)
        {
          // write displacement into global vector
          for (int dim = 0; dim < (3 - isPseudo2D_); dim++)
          {
            int error = -1234;

            error = cell_adhesion_disp_->SumIntoGlobalValue(dofs[dim], 0, adhesiondisp(dim));
            if (error != 0) dserror("SumIntoMyValue returned err=%d", error);
          }  // dim loop
        }    // if owned
      }      // if matched
    }        // if adhesion_node != NULL

    // increment iterator
    if (it != adh_nod_param_coords_in_backgrd_ele_.end()) it++;

  }  // do for all adhesion nodes

  cell_adhesion_disp_->Norm2(&adhesiondispnorm);

  if (myrank_ == 0)
  {
    std::cout << "###   L2-Norm of Cell Adhesion Disp Vector  : " << std::setprecision(11)
              << adhesiondispnorm << std::endl;
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::ReadRestart(int step)
{
  cellscatra_subproblem_->ReadRestart(step);
  poroscatra_subproblem_->ReadRestart(step);
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::PrepareOutput()
{
  cellscatra_subproblem_->StructureField()->PrepareOutput();
  poroscatra_subproblem_->PrepareOutput();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::Output()
{
  cellscatra_subproblem_->StructureField()->Output();
  cellscatra_subproblem_->ScaTraField()->ScaTraField()->Output();
  poroscatra_subproblem_->Output();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::Update()
{
  cellscatra_subproblem_->StructureField()->Update();
  cellscatra_subproblem_->ScaTraField()->ScaTraField()->Update();
  poroscatra_subproblem_->Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::BuildImmersedDirichMap(
    const Teuchos::RCP<const DRT::Discretization>& dis,
    const Teuchos::RCP<const Epetra_Map>& dirichmap_orig, Teuchos::RCP<Epetra_Map>& dirichmap_new,
    std::map<int, int>& adh_nod_backgrd_ele_mapping)
{
  std::map<int, int>::iterator iter;
  std::vector<int> mydirichrowdofs(0);

  for (iter = adh_nod_backgrd_ele_mapping.begin(); iter != adh_nod_backgrd_ele_mapping.end();
       iter++)
  {
    DRT::Node* node_ptr = dis->gNode(iter->first);
    if (node_ptr == NULL) dserror("Could not get node with id %d", iter->first);

    if (node_ptr->Owner() == myrank_)
    {
      std::vector<int> dofs = dis->Dof(0, node_ptr);

      for (int dim = 0; dim < 3; ++dim)
      {
        // if not already in original dirich map
        if (dirichmap_orig->LID(dofs[dim]) == -1) mydirichrowdofs.push_back(dofs[dim]);
      }
    }  // if proc owns node
  }    // for all nodes in provided map

  int nummydirichvals = mydirichrowdofs.size();
  dirichmap_new =
      Teuchos::rcp(new Epetra_Map(-1, nummydirichvals, &(mydirichrowdofs[0]), 0, dis->Comm()));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::DoImmersedDirichletCond(
    const Teuchos::RCP<Epetra_Vector>& statevector, const Teuchos::RCP<Epetra_Vector>& dirichvals,
    const Teuchos::RCP<const Epetra_Map>& dbcmap)
{
  int mynumvals = dbcmap->NumMyElements();
  double* myvals = dirichvals->Values();

  for (int i = 0; i < mynumvals; ++i)
  {
    int gid = dbcmap->GID(i);

#ifdef DEBUG
    if (mynumvals == 0) dserror("dbcmap empty!");
    int err = -2;
    int lid = dirichvals->Map().LID(gid);
    err = statevector->ReplaceGlobalValue(gid, 0, myvals[lid]);
    if (err == -1)
      dserror("VectorIndex >= NumVectors()");
    else if (err == 1)
      dserror("GlobalRow not associated with calling processor");
    else if (err != -1 and err != 1 and err != 0)
      dserror("Trouble using ReplaceGlobalValue on fluid state vector. ErrorCode = %d", err);
#else
    int lid = dirichvals->Map().LID(gid);
    statevector->ReplaceGlobalValue(gid, 0, myvals[lid]);
#endif
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::ResetImmersedInformation()
{
  if (myrank_ == 0) std::cout << "\n Reset Immersed Information ...\n" << std::endl;

  Teuchos::ParameterList params;
  params.set<int>("action", FLD::reset_immersed_ele);
  params.set<int>("Physical Type", poroscatra_subproblem_->FluidField()->PhysicalType());
  backgroundfluiddis_->Evaluate(params);

  immersed_information_invalid_ = true;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedAdhesionTraction::PrintStepInfo()
{
  if (myrank_ == 0)
  {
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
    std::cout << "### DO CELL-ECM ADHESION INTERACTION STEP ...                                    "
                 "            ###"
              << std::endl;
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
  }
}
