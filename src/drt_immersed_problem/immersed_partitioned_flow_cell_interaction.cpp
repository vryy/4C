/*----------------------------------------------------------------------*/
/*! \file

\brief partitioned immersed cell - interstitial flow interaction algorithm

\level 2

\maintainer Jonas Eichinger
*----------------------------------------------------------------------*/
#include "immersed_partitioned_flow_cell_interaction.H"

#include "../linalg/linalg_utils.H"

#include "../drt_poroelast/poro_scatra_base.H"
#include "../drt_poroelast/poroelast_utils_setup.H"

#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_adapter/ad_scatra_wrapper_cellmigration.H"
#include "../drt_adapter/ad_str_multiphysicswrapper_cellmigration.H"

#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../drt_inpar/inpar_immersed.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_structure/stru_aux.H"

#include "../drt_inpar/inpar_cell.H"

#include "../drt_lib/drt_condition.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_mortar/mortar_calc_utils.H"


IMMERSED::ImmersedPartitionedFlowCellInteraction::ImmersedPartitionedFlowCellInteraction(
    const Teuchos::ParameterList& params, const Epetra_Comm& comm)
    : ImmersedPartitioned(comm),
      artificial_velocity_isvalid_(false),
      boundary_traction_isvalid_(false),
      immersed_info_isvalid_(false)
{
  // get pointer to global problem
  globalproblem_ = DRT::Problem::Instance();
  // safety check
  if (globalproblem_->CellMigrationParams().get<std::string>("FLUID_INTERACTION") != "yes")
    dserror("Parameter FLUID_INTERACTION must be set to 'yes' in ---CELL DYNAMIC section.");

  // important variables for parallel simulations
  myrank_ = comm.MyPID();
  numproc_ = comm.NumProc();

  // get pointer to discretizations
  backgroundfluiddis_ = globalproblem_->GetDis("porofluid");
  backgroundstructuredis_ = globalproblem_->GetDis("structure");
  immerseddis_ = globalproblem_->GetDis("cell");
  scatradis_ = globalproblem_->GetDis("scatra");

  // check whether deformable background mesh is requested
  isALE_ =
      (globalproblem_->ImmersedMethodParams().get<std::string>("DEFORM_BACKGROUND_MESH") == "yes");
  if (isALE_ and myrank_ == 0) std::cout << " Deformable background mesh is used." << std::endl;

  // get coupling variable
  displacementcoupling_ = globalproblem_->CellMigrationParams()
                              .sublist("FLOW INTERACTION MODULE")
                              .get<std::string>("COUPVARIABLE") == "Displacement";
  if (displacementcoupling_ and myrank_ == 0)
    std::cout << " Coupling variable for partitioned FSI scheme :  Displacements " << std::endl;
  else if (!displacementcoupling_ and myrank_ == 0)
    std::cout << " Coupling variable for partitioned FSI scheme :  Force " << std::endl;

  // set switch for interface velocity correction
  correct_boundary_velocities_ = (DRT::INPUT::IntegralValue<int>(
      globalproblem_->ImmersedMethodParams(), "CORRECT_BOUNDARY_VELOCITIES"));

  // print acceleration method
  if (globalproblem_->CellMigrationParams()
          .sublist("FLOW INTERACTION MODULE")
          .get<std::string>("COUPALGO") == "iter_stagg_fixed_rel_param")
  {
    is_relaxation_ = false;
    if (myrank_ == 0) std::cout << "\n Using FIXED relaxation parameter. " << std::endl;
  }
  else if (globalproblem_->CellMigrationParams()
               .sublist("FLOW INTERACTION MODULE")
               .get<std::string>("COUPALGO") == "iter_stagg_AITKEN_rel_param")
  {
    is_relaxation_ = true;
    if (myrank_ == 0) std::cout << "\n Using AITKEN relaxation parameter. " << std::endl;
  }
  else
    dserror("Unknown definition of COUPALGO in FLOW INTERACTION MODULE section for Immersed CFI.");

  // check for unfeasible combination
  if (correct_boundary_velocities_ and displacementcoupling_ and is_relaxation_)
    dserror(
        "Interface velocity correction is not possible with displacement coupled Immersed FSI in "
        "combination with relaxation.");

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

  // get pointer to scatra time integration wrapper for the poro-scatra field
  poroscatra_wrapper_ = params.get<Teuchos::RCP<ADAPTER::AdapterScatraWrapperCellMigration>>(
      "RCPToPoroScatraWrapper");

  // get pointer to cell structure
  Teuchos::RCP<ADAPTER::MultiphysicsStructureWrapperCellMigration> multiphysicswrapper =
      params.get<Teuchos::RCP<ADAPTER::MultiphysicsStructureWrapperCellMigration>>(
          "RCPToCellStructure");
  // safety check
  if (multiphysicswrapper == Teuchos::null)
    dserror("no pointer to MultiphysicsStructureWrapperCellMigration provided");
  // get the fsi specific structure wrapper
  cellstructure_ = multiphysicswrapper->GetFSIStructureWrapperPtr();

  // create instance of poroelast subproblem
  poroscatra_subproblem_ = params.get<Teuchos::RCP<POROELAST::PoroScatraBase>>("RCPToPoroScatra");

  // check object pointers
  if (fluid_SearchTree_ == Teuchos::null) dserror("no pointer to fluid_SearchTree_ provided !");
  if (cell_SearchTree_ == Teuchos::null) dserror("no pointer to cell_SearchTree_ provided !");
  if (currpositions_cell_ == NULL) dserror("no pointer to currpositions_cell_ provided !");
  if (currpositions_ECM_ == NULL) dserror("no pointer to currpositions_ECM_ provided !");
  if (cellstructure_ == Teuchos::null) dserror("no pointer to cellstructure_ provided !");
  if (poroscatra_subproblem_ == Teuchos::null)
    dserror("no pointer to poroscatra_subproblem_ provided !");
  if (poroscatra_wrapper_ == Teuchos::null) dserror("no pointer to poroscatra_wrapper_ provided !");

  // counter for continued (unconverged) steps
  continued_steps_ = 0;

  // decide whether multiple cell bodies or not
  std::vector<DRT::Condition*> conditions;
  immerseddis_->GetCondition("ImmersedSearchbox", conditions);
  if ((int)conditions.size() > 0)
  {
    if (myrank_ == 0)
      std::cout << " MULTI CELL MIGRATION SIMULATION   Number of cells: " << (int)conditions.size()
                << std::endl;
    multicellmigration_ = true;
  }
  else
    multicellmigration_ = false;

  // 0 undefined , 1 ameboid , 2 proteolytic
  migrationtype_ =
      DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(), "MIGRATIONTYPE");
  if (migrationtype_ == INPAR::CELL::cell_migration_ameboid and myrank_ == 0)
    std::cout << " AMEBOID TYPE MIGRATION. No proteolytic reaction in ECM." << std::endl;
  if (migrationtype_ == INPAR::CELL::cell_migration_proteolytic and myrank_ == 0)
    std::cout << "MESENCHYMAL TYPE MIGRATION. Proteolytic reaction in ECM." << std::endl;
  else if (migrationtype_ == INPAR::CELL::cell_migration_undefined)
    dserror(
        "set MIGRATIONTYPE to 'ameboid' or 'proteolytic' in --CELL DYNAMICS section in your .dat "
        "file.");

  // initialize segregation variables
  segregationconstant_ = globalproblem_->CellMigrationParams()
                             .sublist("PROTEOLYSIS MODULE")
                             .get<double>("SEGREGATION_CONST");
  segregationtype_ = DRT::INPUT::IntegralValue<int>(
      globalproblem_->CellMigrationParams().sublist("PROTEOLYSIS MODULE"), "SEGREGATION");
  segregationby_ = DRT::INPUT::IntegralValue<int>(
      globalproblem_->CellMigrationParams().sublist("PROTEOLYSIS MODULE"), "SEGREGATION_BY");

  // set pointer to poro fpsi structure
  porostructure_ = poroscatra_subproblem_->PoroField()->StructureField();

  // vector of fluid stresses interpolated to cell bdry int points integrated over cell surface
  cell_bdry_traction_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()), true));
  // vector of pressure force interpolated to cell bdry int points integrated over cell surface
  cell_bdry_traction_pressure_part_ =
      Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()), true));
  // vector with fluid velocities interpolated from structure
  porofluid_artificial_velocity_ =
      Teuchos::rcp(new Epetra_Vector(*(poroscatra_subproblem_->FluidField()->DofRowMap()), true));
  // vector with fluid velocities interpolated from structure
  poroscatra_segregated_phi_ =
      Teuchos::rcp(new Epetra_Vector(*(poroscatra_subproblem_->ScaTraField()->DofRowMap(0)), true));

  // construct immersed exchange manager.
  // singleton class that makes immersed variables accessible from everywhere in the code.
  // IsFluidInteraction() is asked for in fluid_ele_calc_poro_p1_immersed.
  exchange_manager_ = DRT::ImmersedFieldExchangeManager::Instance();
  exchange_manager_->SetIsFluidInteraction(true);
  exchange_manager_->SetIsInitialized(true);

  // PSEUDO2D switch
  isPseudo2D_ = DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(), "PSEUDO2D");

  // get integration rule for fluid elements cut by structural boundary
  int num_gp_fluid_bound = globalproblem_->ImmersedMethodParams().get<int>("NUM_GP_FLUID_BOUND");
  if (num_gp_fluid_bound == 8)
    degree_gp_fluid_bound_ = 3;
  else if (num_gp_fluid_bound == 64)
    degree_gp_fluid_bound_ = 7;
  else if (num_gp_fluid_bound == 125)
    degree_gp_fluid_bound_ = 9;
  else if (num_gp_fluid_bound == 343)
    degree_gp_fluid_bound_ = 13;
  else if (num_gp_fluid_bound == 729)
    degree_gp_fluid_bound_ = 17;
  else if (num_gp_fluid_bound == 1000)
    degree_gp_fluid_bound_ = 19;
  else
    dserror(
        "Invalid value for parameter NUM_GP_FLUID_BOUND (valid parameters are 8, 64, 125, 343, 729 "
        "and 1000). \n"
        "Fix your Input File.");

  // Validation flag for velocity in artificial domain. After each structure solve the velocity
  // becomes invalid and needs to be projected again.
  artificial_velocity_isvalid_ = false;

  // Validation flag for bdry. traction on structure. After each fluid solve the traction becomes
  // invalid and needs to be integrated again.
  boundary_traction_isvalid_ = false;

  // Validation flag for immersed information. After each structure solve, we have to assess
  // again, which fluid elements are covered by the structure and which fluid eles are cut by
  // the interface.
  // todo: NOTE: There is little inconsistency right now in this method.
  //             Fluid elements are labeled isboundarimmersed_ if:
  //             1) at least one but not all nodes are covered, or
  //             2) if a gp of the structural surface lies in a fluid element.
  //
  //             The first criterion is checked in CalcArtificialVelocity(). This is
  //             done at the same time as struct vel. projection is performed.
  //             The second criterion is checked in CalcFluidTractionsOnStructure().
  //             This is done at the same time as the bdry. traction is integrated.
  //
  //             Since, the fluid field decides which gps of the fluid are compressible,
  //             based on IsImmersed() and IsBoundaryImmersed() information, it might
  //             happen, that after performing CalcArtificialVelocity(), the nodal criterion
  //             is updated correctly, but since the structure has moved, the gp criterion
  //             is invalid. It is only updated after the next struct solve. So the fluid
  //             might be solved with incorrect compressible gps.
  //
  //             To fix this, one would have to split the information update and the projections.
  //
  //             However:
  //             1) This would be more expensive.
  //             2) I assume the error we make by the procedure explained above is very small, since
  //                the deformations between iteration steps should be small. Especially, when we
  //                are close to convergence, the difference in structure deformation should be so
  //                small, that no new gps should have changed their states inbetween the last two
  //                structure solves.
  immersed_info_isvalid_ = false;

  // wait for all processors to arrive here
  Comm().Barrier();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int IMMERSED::ImmersedPartitionedFlowCellInteraction::Init(const Teuchos::ParameterList& params)
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
void IMMERSED::ImmersedPartitionedFlowCellInteraction::Setup()
{
  // make sure Init(...) was called first
  CheckIsInit();

  // get parameters for nox
  const Teuchos::ParameterList& noxparams =
      globalproblem_->CellMigrationParams().sublist("FLOW INTERACTION MODULE");
  SetDefaultParameters(noxparams, NOXParameterList());
  // noxparameterlist_.print();

  // set flag issetup true
  SetIsSetup(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::CouplingOp(
    const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  CheckIsInit();
  CheckIsSetup();

  // DISPLACEMENT COUPLING
  if (displacementcoupling_)
  {
    // get the current artificial velocity state
    const Teuchos::RCP<Epetra_Vector> artificial_velocity_n = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL BackgroundOp
    ////////////////////
    BackgroundOp(artificial_velocity_n, fillFlag);  //!< solve the interstitial flow

    ////////////////////
    // CALL ImmersedOp
    ////////////////////
    CalcFluidTractionOnStructure();  //!< // calculate new fluid traction interpolated to structural
                                     //!< surface
    ImmersedOp(cellstructure_->Interface()->ExtractIMMERSEDCondVector(cell_bdry_traction_),
        fillFlag);               //!< solve the cell
    ResetImmersedInformation();  //!< structure moved; immersed info are invalid -> reset
    const Teuchos::RCP<Epetra_Vector> artificial_velocity_np =
        CalcArtificialVelocity();  //!< calc new projected velocities and update immersed
                                   //!< information

    // calc residual
    int err = F.Update(1.0, *artificial_velocity_np, -1.0, *artificial_velocity_n, 0.0);
    if (err != 0) dserror("Vector update of Coupling-residual returned err=%d", err);
  }
  // FORCE COUPLING
  else if (!displacementcoupling_)
  {
    // get the current interface force state
    const Teuchos::RCP<Epetra_Vector> iforcen = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL ImmersedOp
    ////////////////////
    ImmersedOp(iforcen, fillFlag);  //!< solve the cell
    ResetImmersedInformation();     //!< structure moved; immersed info are invalid -> reset

    ////////////////////
    // CALL BackgroundOp
    ////////////////////
    CalcArtificialVelocity();  //!< calc the new velocity in the artificial fluid domain, immersed
                               //!< info are set inside
    BackgroundOp(porofluid_artificial_velocity_, fillFlag);  //!< solve the fluid
    CalcFluidTractionOnStructure();  //!< calculate new fluid traction integrated over structural
                                     //!< surface

    // calc residual
    int err = F.Update(1.0,
        *(cellstructure_->Interface()->ExtractIMMERSEDCondVector(cell_bdry_traction_)), -1.0,
        *iforcen, 0.0);
    if (err != 0) dserror("Vector update of FSI-residual returned err=%d", err);
  }

  // perform n steps max; then set converged
  static bool nlnsolver_continue =
      globalproblem_->ImmersedMethodParams().get<std::string>("DIVERCONT") == "continue";
  static int itemax =
      globalproblem_->CellMigrationParams().sublist("FLOW INTERACTION MODULE").get<int>("ITEMAX");
  if ((IterationCounter())[0] == itemax and nlnsolver_continue)
  {
    double Fnorm = -1234.0;
    F.Norm2(&Fnorm);
    Fnorm = Fnorm / sqrt(F.GlobalLength());

    if (myrank_ == 0)
    {
      std::cout << "\n  Continue with next time step after ITEMAX = " << (IterationCounter())[0]
                << " iterations. L2-Norm of F=" << std::setprecision(12) << Fnorm << "\n"
                << std::endl;
    }
    // set F to zero to tell NOX that this timestep is converged
    F.Scale(0.0);

    // increment the counter for unconverged steps
    continued_steps_++;
  }

  if (globalproblem_->ImmersedMethodParams().get<std::string>("TIMESTATS") == "everyiter")
  {
    Teuchos::TimeMonitor::summarize();
    Teuchos::TimeMonitor::zeroOutTimers();
  }

  // clear states after time step was set converged
  immerseddis_->ClearState();
  backgroundfluiddis_->ClearState();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::BackgroundOp(
    Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values, const FillType fillFlag)
{
  // print
  IMMERSED::ImmersedPartitioned::BackgroundOp(backgrd_dirichlet_values, fillFlag);

  if (fillFlag == User)
  {
    dserror("fillFlag == User : not yet implemented");
  }
  else
  {
    // apply the given artificial velocity to the fluid field
    ApplyImmersedDirichlet(backgrd_dirichlet_values);

    // apply the given artificial velocity to the fluid field
    ApplyImmersedDirichletScatra(poroscatra_segregated_phi_);

    // solve poro
    poroscatra_subproblem_->Solve();

    // remove immersed dirichlets from dbcmap of fluid (may be different in next iteration)
    // RemoveDirichCond();

    // correct the quality of the interface solution
    CorrectInterfaceVelocity();

    // remove scatra values from dbcmap of fluid (may be different in next iteration)
    RemoveDirichCondScatra();

  }  // fillflag is not User

  // we just invalidated the boundary tractions
  boundary_traction_isvalid_ = false;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedFlowCellInteraction::ImmersedOp(
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

    // structure moved; we just invalidated the artificial velocity
    artificial_velocity_isvalid_ = false;

    // we also invalidated the immersed info
    immersed_info_isvalid_ = false;

    return cellstructure_->ExtractImmersedInterfaceDispnp();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  if (myrank_ == 0) std::cout << "Cell Predictor: " << std::endl;
  cellstructure_->PrepareTimeStep();
  if (myrank_ == 0) std::cout << "Poro Predictor: " << std::endl;
  poroscatra_subproblem_->PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedFlowCellInteraction::InitialGuess()
{
  if (displacementcoupling_)
    return CalcArtificialVelocity();
  else  // FORCE COUPLING
    return cellstructure_->Interface()->ExtractIMMERSEDCondVector(cell_bdry_traction_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::BuildImmersedDirichMap(
    Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<Epetra_Map>& dirichmap,
    const Teuchos::RCP<const Epetra_Map>& dirichmap_original, int dofsetnum)
{
  const Epetra_Map* elerowmap = dis->ElementRowMap();
  std::vector<int> mydirichdofs(0);

  for (int i = 0; i < elerowmap->NumMyElements(); ++i)
  {
    // dynamic_cast necessary because virtual inheritance needs runtime information
    DRT::ELEMENTS::FluidImmersedBase* immersedele =
        dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(dis->gElement(elerowmap->GID(i)));
    if (immersedele->HasProjectedDirichlet())
    {
      DRT::Node** nodes = immersedele->Nodes();
      for (int inode = 0; inode < (immersedele->NumNode()); inode++)
      {
        if (static_cast<IMMERSED::ImmersedNode*>(nodes[inode])->IsMatched() and
            nodes[inode]->Owner() == myrank_)
        {
          std::vector<int> dofs = dis->Dof(dofsetnum, nodes[inode]);

          for (int dim = 0; dim < 3; ++dim)
          {
            if (dirichmap_original->LID(dofs[dim]) == -1)  // if not already in original dirich map
              mydirichdofs.push_back(dofs[dim]);
          }
        }
      }
    }
  }

  int nummydirichvals = mydirichdofs.size();
  dirichmap = Teuchos::rcp(new Epetra_Map(-1, nummydirichvals, &(mydirichdofs[0]), 0, dis->Comm()));

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::BuildImmersedScaTraDirichMap(
    Teuchos::RCP<DRT::Discretization> immersedinfodis, Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<Epetra_Map>& dirichmap, const Teuchos::RCP<const Epetra_Map>& dirichmap_original,
    int dofsetnum)
{
  const Epetra_Map* elerowmap = dis->ElementRowMap();
  std::vector<int> mydirichdofs(0);

  for (int i = 0; i < elerowmap->NumMyElements(); ++i)
  {
    // dynamic_cast necessary because virtual inheritance needs runtime information
    int gid = elerowmap->GID(i);
    DRT::ELEMENTS::FluidImmersedBase* immersedele =
        dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(immersedinfodis->gElement(gid));
    if (immersedele->HasProjectedDirichlet())
    {
      DRT::Node** immersednodes = immersedele->Nodes();
      DRT::Node** nodes = dis->gElement(gid)->Nodes();
      for (int inode = 0; inode < (immersedele->NumNode()); inode++)
      {
        if (static_cast<IMMERSED::ImmersedNode*>(immersednodes[inode])->IsMatched() and
            immersednodes[inode]->Owner() == myrank_)
        {
          std::vector<int> dofs = dis->Dof(dofsetnum, nodes[inode]);

          if (dirichmap_original->LID(dofs[0]) == -1)  // if not already in original dirich map
            mydirichdofs.push_back(dofs[0]);
        }
      }
    }
  }

  int nummydirichvals = mydirichdofs.size();
  dirichmap = Teuchos::rcp(new Epetra_Map(-1, nummydirichvals, &(mydirichdofs[0]), 0, dis->Comm()));

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::DoImmersedDirichletCond(
    Teuchos::RCP<Epetra_Vector> statevector, Teuchos::RCP<Epetra_Vector> dirichvals,
    Teuchos::RCP<Epetra_Map> dbcmap)
{
  int mynumvals = dbcmap->NumMyElements();
  double* myvals = dirichvals->Values();

  for (int i = 0; i < mynumvals; ++i)
  {
    int gid = dbcmap->GID(i);

#ifdef DEBUG
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
void IMMERSED::ImmersedPartitionedFlowCellInteraction::PrepareBackgroundOp()
{
  TEUCHOS_FUNC_TIME_MONITOR("IMMERSED::PrepareBackgroundOp()");

  // search radius factor around center of structure bounding box (fac*diagonal of bounding box)
  double structsearchradiusfac =
      DRT::Problem::Instance()->ImmersedMethodParams().get<double>("STRCT_SRCHRADIUS_FAC");

  //  immerseddis_->SetState(0,"displacement",cellstructure_->Dispnp());
  //  immerseddis_->SetState(0,"velocity",cellstructure_->Velnp());
  //  backgroundstructuredis_->SetState(0,"displacement",porostructure_->Dispnp());
  //  backgroundfluiddis_->SetState(0,"veln",poroscatra_subproblem_->FluidField()->Veln());

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
  int procs[numproc_];
  for (int i = 0; i < numproc_; i++) procs[i] = i;
  LINALG::Gather<int, LINALG::Matrix<3, 1>>(
      my_currpositions_cell, *currpositions_cell_, numproc_, &procs[0], Comm());

  // take special care in case of multi-cell migration
  if (multicellmigration_ == false)
  {
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

    if (isALE_) UpdateCurrentPositionsBackgroundNodes();

    SearchPotentiallyCoveredBackgrdElements(&curr_subset_of_backgrounddis_, fluid_SearchTree_,
        *backgroundfluiddis_, *currpositions_ECM_, boundingboxcenter,
        structsearchradiusfac * max_radius, 0);

    if (curr_subset_of_backgrounddis_.empty() == false)
      std::cout << "\nPrepareBackgroundOp returns "
                << curr_subset_of_backgrounddis_.begin()->second.size()
                << " background elements on Proc " << Comm().MyPID() << std::endl;
  }
  else
  {
    dserror("multi cell migration not implemented yet ... ");
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::PrepareOutput()
{
  poroscatra_subproblem_->PrepareOutput();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::Output()
{
  cellstructure_->Output();
  poroscatra_subproblem_->Output();

  WriteDrag();

  if (myrank_ == 0)
    std::cout << " Number of unconverged steps: " << continued_steps_ << "\n" << std::endl;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::WriteDrag()
{
  std::vector<double> drag_total =
      CalcGlobalResultantfromEpetraVector(Comm(), immerseddis_, cell_bdry_traction_);
  std::vector<double> drag_pressure =
      CalcGlobalResultantfromEpetraVector(Comm(), immerseddis_, cell_bdry_traction_pressure_part_);
  std::vector<double> reaction_force =
      CalcGlobalResultantfromEpetraVector(Comm(), immerseddis_, cellstructure_->Freact());

  WriteExtraOutput(Comm(), Time(), "drag", drag_total, drag_pressure, reaction_force);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::Update()
{
  cellstructure_->Update();
  poroscatra_subproblem_->Update();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::ReadRestart(int step)
{
  cellstructure_->ReadRestart(step);
  poroscatra_subproblem_->ReadRestart(step);
  SetTimeStep(poroscatra_subproblem_->PoroField()->Time(), step);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::RemoveDirichCond()
{
  if (dbcmap_immersed_ != Teuchos::null)
    poroscatra_subproblem_->FluidField()->RemoveDirichCond(dbcmap_immersed_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::RemoveDirichCondScatra()
{
  // remove immersed dirichlets from dbcmap of scatra (may be different in next iteration)
  if (migrationtype_ == INPAR::CELL::cell_migration_proteolytic and
      segregationby_ == INPAR::CELL::segregation_by_dirichlet)
  {
    poroscatra_subproblem_->ScaTraField()->RemoveDirichCond(dbcmap_immersed_scatra_);
    // rebuild the combined dbcmap
    poroscatra_subproblem_->BuildCombinedDBCMap();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::CalcFluidTractionOnStructure()
{
  // sanity check
  if (boundary_traction_isvalid_)
    dserror(
        "Boundary traction from fluid onto immersed cell is still valid!\n"
        "If you really need to calc them anew, invalidate flag boundary_traction_isvalid_ at the "
        "proper position.");

  // reinitialize the transfer vector
  cell_bdry_traction_->Scale(0.0);
  cell_bdry_traction_pressure_part_->Scale(0.0);

  // declare and fill parameter list
  Teuchos::ParameterList params;
  params.set<std::string>("action", "calc_fluid_traction");
  params.set<std::string>("backgrddisname", "porofluid");
  params.set<std::string>("immerseddisname", "cell");

  // set the states needed for evaluation
  SetStatesImmersedOP();

  DRT::AssembleStrategy cell_fld_bdry_strategy(0,  // struct dofset for row
      0,                                           // struct dofset for column
      Teuchos::null,                               // matrix 1
      Teuchos::null,                               //
      cell_bdry_traction_,                         // vector 1
      cell_bdry_traction_pressure_part_,           //
      Teuchos::null                                //
  );

  if (myrank_ == 0)
  {
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
    std::cout << "###   Interpolate fluid stresses to structural surface and calculate traction    "
                 "              "
              << std::endl;
  }
  double normofstructbdrytraction = -1234.0;
  EvaluateInterpolationCondition(
      immerseddis_, params, cell_fld_bdry_strategy, "IMMERSEDCoupling", -1);

  cell_bdry_traction_->Norm2(&normofstructbdrytraction);
  if (myrank_ == 0)
  {
    std::cout << "###   Norm of Boundary Fluid Traction:   " << std::setprecision(10)
              << normofstructbdrytraction << std::endl;
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
  }

  // we just validate the boundary tractions
  boundary_traction_isvalid_ = true;

  // we just validated the immersed info again.
  // technically this is not entirely true.
  // see remark in constructor of this class.
  // here we additionally validated the IsBoundaryImmersed
  // information based on the struct. bdry. int. points.
  immersed_info_isvalid_ = true;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFlowCellInteraction::CalcArtificialVelocity()
{
  if (not artificial_velocity_isvalid_)
  {
    // reinitialize the transfer vector
    porofluid_artificial_velocity_->Scale(0.0);

    // remove previously applied immersed Dirichlet conditions
    RemoveDirichCond();

    // declare parameter list
    Teuchos::ParameterList params;

    // provide number of integration points in fluid elements cut by boundary
    params.set<int>("intpoints_fluid_bound", degree_gp_fluid_bound_);
    // provide name of immersed discretization
    params.set<std::string>("immerseddisname", "cell");

    // set the states needed for evaluation
    SetStatesBackgroundOP();
    // update search trees, etc. ...
    PrepareBackgroundOp();

    // calc the fluid velocity from the structural displacements
    DRT::AssembleStrategy fluid_vol_strategy(0,  // struct dofset for row
        0,                                       // struct dofset for column
        Teuchos::null,                           // matrix 1
        Teuchos::null,                           //
        porofluid_artificial_velocity_,          // vector 1
        Teuchos::null,                           //
        Teuchos::null                            //
    );

    if (myrank_ == 0)
    {
      std::cout << "###############################################################################"
                   "#################"
                << std::endl;
      std::cout << "###   Interpolate Dirichlet Values from immersed elements which overlap the "
                << backgroundfluiddis_->Name() << " nodes ..." << std::endl;
      std::cout << "###############################################################################"
                   "#################"
                << std::endl;
    }

    EvaluateImmersed(params, backgroundfluiddis_, &fluid_vol_strategy,
        &curr_subset_of_backgrounddis_, cell_SearchTree_, currpositions_cell_,
        (int)FLD::interpolate_velocity_to_given_point_immersed, false);

    // we just validated the artificial velocity
    artificial_velocity_isvalid_ = true;

    // we just validated the immersed info again.
    // technically this is not entirely true.
    // see remark in constructor of this class.
    immersed_info_isvalid_ = true;
  }

  return porofluid_artificial_velocity_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::ApplyImmersedDirichlet(
    Teuchos::RCP<Epetra_Vector> artificial_velocity)
{
  BuildImmersedDirichMap(backgroundfluiddis_, dbcmap_immersed_,
      poroscatra_subproblem_->FluidField()->GetDBCMapExtractor()->CondMap(), 0);
  poroscatra_subproblem_->FluidField()->AddDirichCond(dbcmap_immersed_);

  // apply immersed dirichlets
  DoImmersedDirichletCond(poroscatra_subproblem_->FluidField()->WriteAccessVelnp(),
      artificial_velocity, dbcmap_immersed_);
  double normofvelocities = -1234.0;
  poroscatra_subproblem_->FluidField()
      ->ExtractVelocityPart(artificial_velocity)
      ->Norm2(&normofvelocities);

  if (myrank_ == 0)
  {
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
    std::cout << "###   Norm of Dirichlet Values:   " << std::setprecision(7) << normofvelocities
              << std::endl;
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::ApplyImmersedDirichletScatra(
    Teuchos::RCP<Epetra_Vector> scatra_values)
{
  if (migrationtype_ == INPAR::CELL::cell_migration_proteolytic)
  {
    if (segregationby_ == INPAR::CELL::segregation_by_dirichlet and
        segregationtype_ == INPAR::CELL::segregation_volumetric)
    {
      // only constant segregation so far
      scatra_values->PutScalar(segregationconstant_);
      BuildImmersedScaTraDirichMap(backgroundfluiddis_, scatradis_, dbcmap_immersed_scatra_,
          poroscatra_subproblem_->ScaTraField()->DirichMaps()->CondMap(), 0);
      poroscatra_subproblem_->ScaTraField()->AddDirichCond(dbcmap_immersed_scatra_);
      DoImmersedDirichletCond(
          poroscatra_subproblem_->ScaTraField()->Phinp(), scatra_values, dbcmap_immersed_scatra_);

      // rebuild the combined dbcmap
      poroscatra_subproblem_->BuildCombinedDBCMap();
    }
    else if (segregationby_ == INPAR::CELL::segregation_by_neumann)
    {
      bool evaluateonlyboundary = false;
      if (segregationtype_ == INPAR::CELL::segregation_volumetric)
        evaluateonlyboundary = false;
      else if (segregationtype_ == INPAR::CELL::segregation_surface)
        evaluateonlyboundary = true;
      else
      {
        if (myrank_ == 0)
          std::cout << "WARNING! Undefined SEGREGATION type! "
                       "Volumetric Segregation is assumed by default."
                    << std::endl;
      }

      Teuchos::ParameterList sparams_struct;

      // provide element parameter list with number of dofset associated with displacement dofs on
      // scatra discretization
      sparams_struct.set<int>("ndsdisp", poroscatra_subproblem_->ScaTraField()->NdsDisp());

      sparams_struct.set<int>("action", (int)SCATRA::calc_immersed_element_source);
      sparams_struct.set<double>("segregation_constant", segregationconstant_);

      // calc the fluid velocity from the cell displacements
      DRT::AssembleStrategy scatra_vol_strategy(0,  // struct dofset for row
          0,                                        // struct dofset for column
          Teuchos::null,                            // matrix 1
          Teuchos::null,                            //
          scatra_values,                            // vector 1
          Teuchos::null,                            //
          Teuchos::null                             //
      );
      if (curr_subset_of_backgrounddis_.empty() == false)
        EvaluateScaTraWithInternalCommunication(scatradis_, backgroundfluiddis_,
            &scatra_vol_strategy, &curr_subset_of_backgrounddis_, cell_SearchTree_,
            currpositions_cell_, sparams_struct, evaluateonlyboundary);

      poroscatra_wrapper_->AddContributionToRHS(scatra_values);
    }
    else
    {
      if (migrationtype_ == INPAR::CELL::cell_migration_proteolytic)
        dserror("combination of SEGREGATION parameters is not implemented. Fix your input file.");
    }

    double normofconcentrations = -1234.0;
    scatra_values->Norm2(&normofconcentrations);

    if (myrank_ == 0)
    {
      std::cout << "###############################################################################"
                   "#################"
                << std::endl;
      std::cout << "###   Norm of Segregated Species Dirichlet values:   " << std::setprecision(7)
                << normofconcentrations << std::endl;
      std::cout << "###############################################################################"
                   "#################"
                << std::endl;
    }
  }  // if proteolytic
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::CorrectInterfaceVelocity()
{
  //***********************************************************************************
  // Correct velocity at nodes of fluid elements being cut by the structural surface
  // (fluid is solved a second time with different Dirichlet values)
  //***********************************************************************************
  if (correct_boundary_velocities_)
  {
    // declare parameter list
    Teuchos::ParameterList params;

    // provide number of integration points in fluid elements cut by boundary
    params.set<int>("intpoints_fluid_bound", degree_gp_fluid_bound_);

    // calc the fluid velocity from the structural displacements
    DRT::AssembleStrategy fluid_vol_strategy(0,  // struct dofset for row
        0,                                       // struct dofset for column
        Teuchos::null,                           // matrix 1
        Teuchos::null,                           //
        porofluid_artificial_velocity_,          // vector 1
        Teuchos::null,                           //
        Teuchos::null                            //
    );

    SetStatesVelocityCorrection();

    if (myrank_ == 0)
    {
      std::cout << "\nCorrection step " << std::endl;
      std::cout << "###############################################################################"
                   "#################"
                << std::endl;
      std::cout << "###   Correct Velocity in fluid boundary elements " << std::endl;
    }

    // calculate new dirichlet velocities for fluid elements cut by structure
    EvaluateImmersed(params, backgroundfluiddis_, &fluid_vol_strategy,
        &curr_subset_of_backgrounddis_, cell_SearchTree_, currpositions_cell_,
        (int)FLD::correct_immersed_fluid_bound_vel, true);

    // apply corrected dirichlet values
    DoImmersedDirichletCond(poroscatra_subproblem_->FluidField()->WriteAccessVelnp(),
        porofluid_artificial_velocity_, dbcmap_immersed_);
    double normofvelocities = -1234.0;
    poroscatra_subproblem_->FluidField()
        ->ExtractVelocityPart(porofluid_artificial_velocity_)
        ->Norm2(&normofvelocities);

    if (myrank_ == 0)
    {
      std::cout << "###############################################################################"
                   "#################"
                << std::endl;
      std::cout << "###   Norm of Dirichlet Values:   " << std::setprecision(7) << normofvelocities
                << std::endl;
      std::cout << "###############################################################################"
                   "#################"
                << std::endl;
    }

    // solve fluid again with new dirichlet values
    // solve poro
    poroscatra_subproblem_->Solve();

  }  // correct_boundary_velocities_ finished

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::ResetImmersedInformation()
{
  if (immersed_info_isvalid_)
    dserror(
        "Immersed information are valid! Reconsider your call to ResetImmersedInformation().\n"
        "Did you forget to invalidate the flag immersed_info_isvalid_?");

  if (myrank_ == 0) std::cout << "\nReset Immersed Information ...\n" << std::endl;

  // reset element and node information about immersed method
  Teuchos::ParameterList params;
  params.set<int>("action", FLD::reset_immersed_ele);
  params.set<int>("Physical Type", poroscatra_subproblem_->FluidField()->PhysicalType());
  params.set<int>("intpoints_fluid_bound", degree_gp_fluid_bound_);
  EvaluateSubsetElements(
      params, backgroundfluiddis_, curr_subset_of_backgrounddis_, (int)FLD::reset_immersed_ele);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::SetStatesBackgroundOP()
{
  immerseddis_->SetState(0, "displacement", cellstructure_->Dispnp());
  immerseddis_->SetState(0, "velocity", cellstructure_->Velnp());

  if (isALE_)
    backgroundfluiddis_->SetState(0, "dispnp", poroscatra_subproblem_->FluidField()->Dispnp());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::SetStatesVelocityCorrection()
{
  immerseddis_->SetState(0, "displacement", cellstructure_->Dispnp());
  immerseddis_->SetState(0, "velocity", cellstructure_->Velnp());
  backgroundfluiddis_->SetState(0, "velnp", poroscatra_subproblem_->FluidField()->Velnp());

  if (isALE_)
    backgroundfluiddis_->SetState(0, "dispnp", poroscatra_subproblem_->FluidField()->Dispnp());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::SetStatesImmersedOP()
{
  backgroundfluiddis_->SetState(0, "velnp", poroscatra_subproblem_->FluidField()->Velnp());
  backgroundfluiddis_->SetState(0, "veln", poroscatra_subproblem_->FluidField()->Veln());
  immerseddis_->SetState(0, "displacement", cellstructure_->Dispnp());

  if (isALE_)
    backgroundfluiddis_->SetState(0, "dispnp", poroscatra_subproblem_->FluidField()->Dispnp());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::UpdateCurrentPositionsBackgroundNodes()
{
  // get displacement state
  // set and get displacement state (export happens inside)
  backgroundfluiddis_->SetState("dispnp", poroscatra_subproblem_->FluidField()->Dispnp());
  Teuchos::RCP<const Epetra_Vector> displacements = backgroundfluiddis_->GetState("dispnp");

  int nummyrownodes = backgroundfluiddis_->NumMyRowNodes();

  // update positions
  for (int lid = 0; lid < nummyrownodes; ++lid)
  {
    const DRT::Node* node = backgroundfluiddis_->lRowNode(lid);
    LINALG::Matrix<3, 1> currpos;
    std::vector<int> dofstoextract(4);
    std::vector<double> mydisp(4);

    // get the current displacement
    backgroundfluiddis_->Dof(0, node, 0, dofstoextract);
    DRT::UTILS::ExtractMyValues(*displacements, mydisp, dofstoextract);

    currpos(0) = node->X()[0] + mydisp.at(0);
    currpos(1) = node->X()[1] + mydisp.at(1);
    currpos(2) = node->X()[2] + mydisp.at(2);

    (*currpositions_ECM_)[node->Id()] = currpos;
  }

  return;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitionedFlowCellInteraction::ReturnCouplingInfo()
{
  return cell_bdry_traction_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::PrintStepInfo()
{
  if (myrank_ == 0)
  {
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
    std::cout << "### DO CELL-FLOW INTERACTION STEP ...                                            "
                 "            ###"
              << std::endl;
    std::cout << "#################################################################################"
                 "###############"
              << std::endl;
  }
}
