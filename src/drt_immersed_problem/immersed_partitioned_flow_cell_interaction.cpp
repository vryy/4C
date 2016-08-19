/*!----------------------------------------------------------------------
\file immersed_partitioned_flow_cell_interaction.cpp

\brief partitioned immersed cell - interstitial flow interaction algorithm

\level 2

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

*----------------------------------------------------------------------*/
#include "immersed_partitioned_flow_cell_interaction.H"

#include "../linalg/linalg_utils.H"

#include "../drt_poroelast/poro_scatra_base.H"
#include "../drt_poroelast/poroelast_utils_setup.H"

#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"

#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../drt_inpar/inpar_immersed.H"

#include "../drt_fsi/fsi_nox_aitken.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_structure/stru_aux.H"

#include "../drt_inpar/inpar_cell.H"

#include "../drt_lib/drt_condition.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_mortar/mortar_calc_utils.H"


IMMERSED::ImmersedPartitionedFlowCellInteraction::ImmersedPartitionedFlowCellInteraction(const Teuchos::ParameterList& params, const Epetra_Comm& comm)
  : ImmersedPartitioned(comm)
{
  // get pointer to fluid search tree from ParameterList
  fluid_SearchTree_ = params.get<Teuchos::RCP<GEO::SearchTree> >("RCPToFluidSearchTree");

  // get pointer to cell search tree from ParameterList
  cell_SearchTree_ = params.get<Teuchos::RCP<GEO::SearchTree> >("RCPToCellSearchTree");

  // get pointer to the current position map of the cell
  currpositions_cell_ = params.get<std::map<int,LINALG::Matrix<3,1> >* >("PointerToCurrentPositionsCell");

  // get pointer to the current position map of the cell
  currpositions_ECM_ = params.get<std::map<int,LINALG::Matrix<3,1> >* >("PointerToCurrentPositionsECM");

  // get pointer to cell structure
  cellstructure_=params.get<Teuchos::RCP<ADAPTER::FSIStructureWrapperImmersed> >("RCPToCellStructure");

  // create instance of poroelast subproblem
  poroscatra_subproblem_ = params.get<Teuchos::RCP<POROELAST::PoroScatraBase> >("RCPToPoroScatra");

  // check object pointers
  if(fluid_SearchTree_==Teuchos::null)
    dserror("no pointer to fluid_SearchTree_ provided !");
  if(cell_SearchTree_==Teuchos::null)
    dserror("no pointer to cell_SearchTree_ provided !");
  if(currpositions_cell_==NULL)
    dserror("no pointer to currpositions_cell_ provided !");
  if(currpositions_ECM_==NULL)
    dserror("no pointer to currpositions_ECM_ provided !");
  if(cellstructure_==Teuchos::null)
    dserror("no pointer to cellstructure_ provided !");
  if(poroscatra_subproblem_==Teuchos::null)
    dserror("no pointer to poroscatra_subproblem_ provided !");

  // important variables for parallel simulations
  myrank_  = comm.MyPID();
  numproc_ = comm.NumProc();

  // get pointer to global problem
  globalproblem_ = DRT::Problem::Instance();

  // safety check
  if (globalproblem_->CellMigrationParams().get<std::string>("FLUID_INTERACTION") != "yes")
    dserror("Parameter FLUID_INTERACTION must be set to 'yes' in ---CELL DYNAMIC section.");

  // get pointer to discretizations
  backgroundfluiddis_     = globalproblem_->GetDis("porofluid");
  backgroundstructuredis_ = globalproblem_->GetDis("structure");
  immerseddis_            = globalproblem_->GetDis("cell");
  scatradis_              = globalproblem_->GetDis("scatra");

  // counter for continued (unconverged) steps
  continued_steps_=0;

  // decide whether multiple cell bodies or not
  std::vector<DRT::Condition*> conditions;
  immerseddis_->GetCondition("ImmersedSearchbox",conditions);
  if((int)conditions.size()>0)
  {
    if(myrank_==0)
      std::cout<<" MULTI CELL MIGRATION SIMULATION   Number of cells: "<<(int)conditions.size()<<std::endl;
    multicellmigration_ = true;
  }
  else
    multicellmigration_ = false;

  // 0 undefined , 1 ameboid , 2 proteolytic
  migrationtype_=DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(),"MIGRATIONTYPE");
  if(migrationtype_==INPAR::CELL::cell_migration_ameboid and myrank_==0)
    std::cout<<"AMEBOID TYPE MIGRATION. No proteolytic reaction in ECM."<<std::endl;
  if(migrationtype_==INPAR::CELL::cell_migration_proteolytic and myrank_==0)
    std::cout<<"MESENCHYMAL TYPE MIGRATION. Proteolytic reaction in ECM."<<std::endl;
  else if(migrationtype_==INPAR::CELL::cell_migration_undefined)
    dserror("set MIGRATIONTYPE to 'ameboid' or 'proteolytic' in --CELL DYNAMICS section in your .dat file.");

  // initialize segregation variables
  segregationconstant_=globalproblem_->CellMigrationParams().sublist("PROTEOLYSIS MODULE").get<double>("SEGREGATION_CONST");
  segregationtype_=DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams().sublist("PROTEOLYSIS MODULE"),"SEGREGATION");
  segregationby_=DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams().sublist("PROTEOLYSIS MODULE"),"SEGREGATION_BY");

  // setup the relaxation parameters
  SetupRelaxation();

  // set pointer to poro fpsi structure
  porostructure_ = poroscatra_subproblem_->PoroField()->StructureField();

  // vector of fluid stresses interpolated to cell bdry int points integrated over cell surface
  cell_bdry_traction_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()),true));
  // vector with fluid velocities interpolated from structure
  porofluid_artificial_velocity_ = Teuchos::rcp(new Epetra_Vector(*(poroscatra_subproblem_->FluidField()->DofRowMap()),true));
  // vector with fluid velocities interpolated from structure
  poroscatra_segregated_phi_ = Teuchos::rcp(new Epetra_Vector(*(poroscatra_subproblem_->ScaTraField()->DofRowMap(0)),true));

  // get coupling variable
  displacementcoupling_ = globalproblem_->ImmersedMethodParams().sublist("PARTITIONED SOLVER").get<std::string>("COUPVARIABLE_FSI") == "Displacement";
  if(displacementcoupling_ and myrank_==0)
    std::cout<<" Coupling variable for partitioned scheme :  Displacements "<<std::endl;
  else if (!displacementcoupling_ and myrank_==0)
    std::cout<<" Coupling variable for partitioned scheme :  Force "<<std::endl;

  // construct immersed exchange manager. singleton class that makes immersed variables comfortably accessible from everywhere in the code
  exchange_manager_ = DRT::ImmersedFieldExchangeManager::Instance();
  exchange_manager_->SetIsFluidInteraction(true);
  exchange_manager_->SetPointerToCurrentSubsetOfBackgrdDis(&curr_subset_of_backgrounddis_);
  exchange_manager_->SetIsInitialized(true);

  // PSEUDO2D switch
  isPseudo2D_ = DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(),"PSEUDO2D");

  // get integration rule for fluid elements cut by structural boundary
  int num_gp_fluid_bound = globalproblem_->ImmersedMethodParams().get<int>("NUM_GP_FLUID_BOUND");
  if(num_gp_fluid_bound == 8)
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
    dserror("Invalid value for parameter NUM_GP_FLUID_BOUND (valid parameters are 8, 64, 125, 343, 729 and 1000). \n"
            "Fix your Input File.");

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::CouplingOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  ReinitTransferVectors();

  ResetImmersedInformation();

  // DISPLACEMENT COUPLING
  if (displacementcoupling_)
  {
    const Teuchos::RCP<Epetra_Vector> idispn = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL BackgroundOp
    ////////////////////
    PrepareBackgroundOp();
    BackgroundOp(porofluid_artificial_velocity_, fillFlag);

    ////////////////////
    // CALL ImmersedOp
    ////////////////////
    PrepareImmersedOp();
    Teuchos::RCP<Epetra_Vector> idispnp =
    ImmersedOp(cellstructure_->Interface()->ExtractIMMERSEDCondVector(cell_bdry_traction_), fillFlag);

    int err = F.Update(1.0, *idispnp, -1.0, *idispn, 0.0);
    if(err != 0)
      dserror("Vector update of Coupling-residual returned err=%d",err);

  }
  // FORCE COUPLING
  else if(!displacementcoupling_)
  {
    const Teuchos::RCP<Epetra_Vector> iforcen = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL ImmersedOp
    ////////////////////
    PrepareImmersedOp();
    ImmersedOp(iforcen, fillFlag);

    ////////////////////
    // CALL BackgroundOp
    ////////////////////
    PrepareBackgroundOp(); // determine elements cut by the boundary
    BackgroundOp(porofluid_artificial_velocity_, fillFlag);

    int err = F.Update(1.0, *(cellstructure_->Interface()->ExtractIMMERSEDCondVector(cell_bdry_traction_)), -1.0, *iforcen, 0.0);

    if(err != 0)
      dserror("Vector update of FSI-residual returned err=%d",err);

  }

  // perform n steps max; then set converged
  static bool nlnsolver_continue = globalproblem_->ImmersedMethodParams().get<std::string>("DIVERCONT") == "continue";
  static int  itemax = globalproblem_->ImmersedMethodParams().sublist("PARTITIONED SOLVER").get<int>("ITEMAX");
  if((IterationCounter())[0] == itemax and nlnsolver_continue)
  {
    if(Comm().MyPID()==0)
      std::cout<<"\n  Continue with next time step after ITEMAX = "<<(IterationCounter())[0]<<" iterations. \n"<<std::endl;

    // !!! EXPERIMENTAL !!!
    // set F to zero to tell NOX that this timestep is converged
    Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp( new Epetra_Vector(F.Map(),true));
    F.Update(1.0,*zeros,0.0);
    // !!! EXPERIMENTAL !!!

    immerseddis_->ClearState();
    backgroundfluiddis_->ClearState();

    continued_steps_++;
  }

  if(DRT::Problem::Instance()->ImmersedMethodParams().get<std::string>("TIMESTATS")=="everyiter")
  {
    Teuchos::TimeMonitor::summarize();
    Teuchos::TimeMonitor::zeroOutTimers();
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::BackgroundOp(Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values,
                                                                    const FillType fillFlag)
{
  // print
  IMMERSED::ImmersedPartitioned::BackgroundOp(backgrd_dirichlet_values,fillFlag);
  Teuchos::ParameterList params;
  params.set<int>("intpoints_fluid_bound",degree_gp_fluid_bound_);

  if (fillFlag==User)
  {
    dserror("fillFlag == User : not yet implemented");
  }
  else
  {
    // calc the fluid velocity from the cell displacements
    DRT::AssembleStrategy fluid_vol_strategy(
        0,              // struct dofset for row
        0,              // struct dofset for column
        Teuchos::null,  // matrix 1
        Teuchos::null,  //
        backgrd_dirichlet_values ,  // vector 1
        Teuchos::null,  //
        Teuchos::null   //
    );

    if(myrank_ == 0)
    {
      std::cout<<"################################################################################################"<<std::endl;
      std::cout<<"###   Interpolate transfer quantities from immersed elements which overlap the "<<backgroundfluiddis_->Name()<<" nodes ..."<<std::endl;
    }

    if(curr_subset_of_backgrounddis_.empty()==false)
      EvaluateImmersed(params,
          backgroundfluiddis_,
          &fluid_vol_strategy,
          &curr_subset_of_backgrounddis_,
          cell_SearchTree_, currpositions_cell_,
          (int)FLD::interpolate_velocity_to_given_point_immersed,
          false);
    else
    {
      // do nothing : a proc without subset of backgrd dis. does not need to enter evaluation,
      //              since cell is ghosted on all procs and no communication has to be performed.
    }

    /////////////////////////////////////////////////////////////////////////////////////
    //  Apply Dirichlet values to background fluid
    //
    BuildImmersedDirichMap(backgroundfluiddis_, dbcmap_immersed_, poroscatra_subproblem_->FluidField()->GetDBCMapExtractor()->CondMap(),0);
    poroscatra_subproblem_->FluidField()->AddDirichCond(dbcmap_immersed_);
    DoImmersedDirichletCond(poroscatra_subproblem_->FluidField()->WriteAccessVelnp(),backgrd_dirichlet_values, dbcmap_immersed_);

    /////////////////////////////////////////////////////////////////////////////////////
    //  Apply Dirichlet values to background scatra field
    //
    if(migrationtype_==INPAR::CELL::cell_migration_proteolytic and segregationby_==INPAR::CELL::segregation_by_dirichlet and segregationtype_==INPAR::CELL::segregation_volumetric)
    {
      // only constant segregation so far
      poroscatra_segregated_phi_->PutScalar(segregationconstant_);
      BuildImmersedScaTraDirichMap(backgroundfluiddis_, scatradis_, dbcmap_immersed_scatra_, poroscatra_subproblem_->ScaTraField()->DirichMaps()->CondMap(),0);
      poroscatra_subproblem_->ScaTraField()->AddDirichCond(dbcmap_immersed_scatra_);
      DoImmersedDirichletCond(poroscatra_subproblem_->ScaTraField()->Phinp(),poroscatra_segregated_phi_, dbcmap_immersed_scatra_);
    }
    /////////////////////////////////////////////////////////////////////////////////////
    //  Apply Source values to rhs of background scatra field
    //
    else if(migrationtype_==INPAR::CELL::cell_migration_proteolytic and segregationby_==INPAR::CELL::segregation_by_neumann)
    {
      bool evaluateonlyboundary = false;
      if(segregationtype_==INPAR::CELL::segregation_volumetric)
        evaluateonlyboundary = false;
      else if(segregationtype_==INPAR::CELL::segregation_surface)
        evaluateonlyboundary = true;
      else
      {
        if(myrank_==0)
          std::cout<<"WARNING! Undefined SEGREGATION type! "
          "Volumetric Segregation is assumed by default."<<std::endl;
      }

      Teuchos::ParameterList sparams_struct;

      // provide element parameter list with number of dofset associated with displacement dofs on scatra discretization
      sparams_struct.set<int>("ndsdisp",poroscatra_subproblem_->ScaTraField()->NdsDisp());

      sparams_struct.set<int>("action",(int)SCATRA::calc_immersed_element_source);
      sparams_struct.set<double>("segregation_constant",segregationconstant_);

      // calc the fluid velocity from the cell displacements
      DRT::AssembleStrategy scatra_vol_strategy(
          0,              // struct dofset for row
          0,              // struct dofset for column
          Teuchos::null,  // matrix 1
          Teuchos::null,  //
          poroscatra_segregated_phi_ ,  // vector 1
          Teuchos::null,  //
          Teuchos::null   //
      );
      if(curr_subset_of_backgrounddis_.empty()==false)
        EvaluateScaTraWithInternalCommunication(scatradis_,
            backgroundfluiddis_,
            &scatra_vol_strategy,
            &curr_subset_of_backgrounddis_,
            cell_SearchTree_,
            currpositions_cell_,
            sparams_struct,
            evaluateonlyboundary);

      poroscatra_subproblem_->ScaTraField()->AddContributionToRHS(poroscatra_segregated_phi_);

    }
    else
    {
      if(migrationtype_==INPAR::CELL::cell_migration_proteolytic)
        dserror("combination of SEGREGATION parameters is not implemented. Fix your input file.");
    }

    // rebuild the combined dbcmap
    poroscatra_subproblem_->BuildCombinedDBCMap();

    double normofvelocities = -1234.0;
    poroscatra_subproblem_->FluidField()->ExtractVelocityPart(backgrd_dirichlet_values)->Norm2(&normofvelocities);

    double normofconcentrations = -1234.0;
    poroscatra_segregated_phi_->Norm2(&normofconcentrations);

    if(myrank_ == 0)
    {
      std::cout<<"###   Norm of Fluid Velocity Dirichlet values:           "<<std::setprecision(7)<<normofvelocities<<std::endl;
      if(migrationtype_==INPAR::CELL::cell_migration_proteolytic)
      {
        std::cout<<"###   Norm of Segregated Species Dirichlet values:     "<<std::setprecision(10)<<normofconcentrations<<std::endl;
      }
      std::cout<<"###   Norm of transferred source/sink:                 "<<std::setprecision(10)<<"not available"<<std::endl;
      std::cout<<"################################################################################################"<<std::endl;
    }

    exchange_manager_->SetImmersedSearchTree(cell_SearchTree_);
    exchange_manager_->SetCurrentPositionsImmersedDis(currpositions_cell_);

    // solve poro
    poroscatra_subproblem_->Solve();

    //std::cout<<"Matched "<<exchange_manager_->GetCheckCounter()<<" integration points in "<<exchange_manager_->GetNumIsImmersedBoundary()<<" IsImmersedBoundary() elements on PROC "<<myrank_<<". "<<std::endl;

    // remove immersed dirichlets from dbcmap of fluid (may be different in next iteration)
    poroscatra_subproblem_->FluidField()->RemoveDirichCond(dbcmap_immersed_);
    // remove immersed dirichlets from dbcmap of scatra (may be different in next iteration)
    if(migrationtype_==INPAR::CELL::cell_migration_proteolytic and segregationby_==INPAR::CELL::segregation_by_dirichlet)
      poroscatra_subproblem_->ScaTraField()->RemoveDirichCond(dbcmap_immersed_scatra_);

    // rebuild the combined dbcmap
    poroscatra_subproblem_->BuildCombinedDBCMap();

  } // fillflag not User


  // calculate new fluid traction interpolated to structural surface
  CalcFluidTractionOnStructure();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFlowCellInteraction::ImmersedOp(Teuchos::RCP<Epetra_Vector> bdry_traction,
                                                       const FillType fillFlag)
{
  IMMERSED::ImmersedPartitioned::ImmersedOp(bdry_traction,fillFlag);

  if(fillFlag==User)
  {
    dserror("fillFlag == User : not yet implemented");
    return Teuchos::null;
  }
  else
  {
    // prescribe neumann values at structural boundary dofs
    cellstructure_->ApplyImmersedInterfaceForcesTemporaryImplementation(Teuchos::null,bdry_traction);

    // solve cell
    cellstructure_->Solve();

    return cellstructure_->ExtractImmersedInterfaceDispnp();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  if(myrank_==0)
    std::cout<<"Cell Predictor: "<<std::endl;
  cellstructure_->PrepareTimeStep();
  if(myrank_==0)
    std::cout<<"Poro Predictor: "<<std::endl;
  poroscatra_subproblem_->PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedFlowCellInteraction::InitialGuess()
{
  if(displacementcoupling_)
    return cellstructure_->PredictImmersedInterfaceDispnp();
  else // FORCE COUPLING
  {
    return cellstructure_->Interface()->ExtractIMMERSEDCondVector(cell_bdry_traction_);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::BuildImmersedDirichMap(Teuchos::RCP<DRT::Discretization> dis,
                                                                        Teuchos::RCP<Epetra_Map>& dirichmap,
                                                                        const Teuchos::RCP<const Epetra_Map>& dirichmap_original,
                                                                        int dofsetnum)
{
  const Epetra_Map* elerowmap = dis->ElementRowMap();
  std::vector<int> mydirichdofs(0);

  for(int i=0; i<elerowmap->NumMyElements(); ++i)
  {
    // dynamic_cast necessary because virtual inheritance needs runtime information
    DRT::ELEMENTS::FluidImmersedBase* immersedele = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(dis->gElement(elerowmap->GID(i)));
    if(immersedele->HasProjectedDirichlet())
    {
      DRT::Node** nodes = immersedele->Nodes();
      for (int inode=0; inode<(immersedele->NumNode()); inode++)
      {
        if(static_cast<IMMERSED::ImmersedNode* >(nodes[inode])->IsMatched() and nodes[inode]->Owner()==myrank_)
        {
          std::vector<int> dofs = dis->Dof(dofsetnum,nodes[inode]);

          for (int dim=0;dim<3;++dim)
          {
            if(dirichmap_original->LID(dofs[dim]) == -1) // if not already in original dirich map
              mydirichdofs.push_back(dofs[dim]);
          }

        }
      }
    }
  }

  int nummydirichvals = mydirichdofs.size();
  dirichmap = Teuchos::rcp( new Epetra_Map(-1,nummydirichvals,&(mydirichdofs[0]),0,dis->Comm()) );

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::BuildImmersedScaTraDirichMap(Teuchos::RCP<DRT::Discretization> immersedinfodis,
                                                                        Teuchos::RCP<DRT::Discretization> dis,
                                                                        Teuchos::RCP<Epetra_Map>& dirichmap,
                                                                        const Teuchos::RCP<const Epetra_Map>& dirichmap_original,
                                                                        int dofsetnum)
{
  const Epetra_Map* elerowmap = dis->ElementRowMap();
  std::vector<int> mydirichdofs(0);

  for(int i=0; i<elerowmap->NumMyElements(); ++i)
  {
    // dynamic_cast necessary because virtual inheritance needs runtime information
    int gid = elerowmap->GID(i);
    DRT::ELEMENTS::FluidImmersedBase* immersedele = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(immersedinfodis->gElement(gid));
    if(immersedele->HasProjectedDirichlet())
    {
      DRT::Node** immersednodes = immersedele->Nodes();
      DRT::Node** nodes = dis->gElement(gid)->Nodes();
      for (int inode=0; inode<(immersedele->NumNode()); inode++)
      {
        if(static_cast<IMMERSED::ImmersedNode* >(immersednodes[inode])->IsMatched() and immersednodes[inode]->Owner()==myrank_)
        {
          std::vector<int> dofs = dis->Dof(dofsetnum,nodes[inode]);

            if(dirichmap_original->LID(dofs[0]) == -1) // if not already in original dirich map
              mydirichdofs.push_back(dofs[0]);

        }
      }
    }
  }

  int nummydirichvals = mydirichdofs.size();
  dirichmap = Teuchos::rcp( new Epetra_Map(-1,nummydirichvals,&(mydirichdofs[0]),0,dis->Comm()) );

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::DoImmersedDirichletCond(Teuchos::RCP<Epetra_Vector> statevector, Teuchos::RCP<Epetra_Vector> dirichvals, Teuchos::RCP<Epetra_Map> dbcmap)
{
  int mynumvals = dbcmap->NumMyElements();
  double* myvals = dirichvals->Values();

  for(int i=0;i<mynumvals;++i)
  {
    int gid = dbcmap->GID(i);

#ifdef DEBUG
    int err = -2;
    int lid = dirichvals->Map().LID(gid);
    err = statevector -> ReplaceGlobalValue(gid,0,myvals[lid]);
    if(err==-1)
      dserror("VectorIndex >= NumVectors()");
    else if (err==1)
        dserror("GlobalRow not associated with calling processor");
    else if (err != -1 and err != 1 and err != 0)
      dserror("Trouble using ReplaceGlobalValue on fluid state vector. ErrorCode = %d",err);
#else
    int lid = dirichvals->Map().LID(gid);
    statevector -> ReplaceGlobalValue(gid,0,myvals[lid]);
#endif

  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::PrepareImmersedOp()
{
  immerseddis_->SetState(0,"velnp",cellstructure_->Velnp());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::PrepareBackgroundOp()
{
  TEUCHOS_FUNC_TIME_MONITOR("IMMERSED::PrepareBackgroundOp()");

  immerseddis_->SetState(0,"displacement",cellstructure_->Dispnp());
  immerseddis_->SetState(0,"velocity",cellstructure_->Velnp());
  backgroundstructuredis_->SetState(0,"displacement",porostructure_->Dispnp());
  backgroundfluiddis_->SetState(0,"veln",poroscatra_subproblem_->FluidField()->Veln());

  double structsearchradiusfac = DRT::Problem::Instance()->ImmersedMethodParams().get<double>("STRCT_SRCHRADIUS_FAC");

  // determine subset of fluid discretization which is potentially underlying the immersed discretization
  //
  // get state
  Teuchos::RCP<const Epetra_Vector> celldisplacements = cellstructure_->Dispnp();

  // find current positions for immersed structural discretization
  std::map<int,LINALG::Matrix<3,1> > my_currpositions_cell;
  for (int lid = 0; lid < immerseddis_->NumMyRowNodes(); ++lid)
  {
    const DRT::Node* node = immerseddis_->lRowNode(lid);
    LINALG::Matrix<3,1> currpos;
    std::vector<int> dofstoextract(3);
    std::vector<double> mydisp(3);

    // get the current displacement
    immerseddis_->Dof(node,0,dofstoextract);
    DRT::UTILS::ExtractMyValues(*celldisplacements,mydisp,dofstoextract);

    currpos(0) = node->X()[0]+mydisp.at(0);
    currpos(1) = node->X()[1]+mydisp.at(1);
    currpos(2) = node->X()[2]+mydisp.at(2);

    my_currpositions_cell[node->Id()] = currpos;
  }

  // communicate local currpositions:
  // map with current cell positions should be same on all procs
  // to make use of the advantages of ghosting the cell redundantly
  // on all procs.
  int procs[Comm().NumProc()];
  for(int i=0;i<Comm().NumProc();i++)
    procs[i]=i;
  LINALG::Gather<int,LINALG::Matrix<3,1> >(my_currpositions_cell,*currpositions_cell_,Comm().NumProc(),&procs[0],Comm());

  if (multicellmigration_ == false)
  {
    // get bounding box of current configuration of structural dis
    const LINALG::Matrix<3,2> structBox = GEO::getXAABBofDis(*immerseddis_,*currpositions_cell_);
    double max_radius = sqrt(pow(structBox(0,0)-structBox(0,1),2)+pow(structBox(1,0)-structBox(1,1),2)+pow(structBox(2,0)-structBox(2,1),2));
    // search for background elements within a certain radius around the center of the immersed bounding box
    LINALG::Matrix<3,1> boundingboxcenter;
    boundingboxcenter(0) = structBox(0,0)+(structBox(0,1)-structBox(0,0))*0.5;
    boundingboxcenter(1) = structBox(1,0)+(structBox(1,1)-structBox(1,0))*0.5;
    boundingboxcenter(2) = structBox(2,0)+(structBox(2,1)-structBox(2,0))*0.5;

    // search background elements covered by bounding box of cell
    curr_subset_of_backgrounddis_ = fluid_SearchTree_->searchElementsInRadius(*backgroundfluiddis_,*currpositions_ECM_,boundingboxcenter,structsearchradiusfac*max_radius,0);

    if(curr_subset_of_backgrounddis_.empty() == false)
      std::cout<<"\nPrepareBackgroundOp returns "<<curr_subset_of_backgrounddis_.begin()->second.size()<<" background elements on Proc "<<Comm().MyPID()<<std::endl;

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

  if(myrank_==0)
    std::cout<<" Number of unconverged steps: "<<continued_steps_<<"\n"<<std::endl;

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
void IMMERSED::ImmersedPartitionedFlowCellInteraction::SetupRelaxation()
{
  // initialize some relaxation related member variables
    relaxforceglobally_ = globalproblem_->ImmersedMethodParams().get<std::string>("APPLY_FORCE_RELAX") == "globally";
    relaxvelglobally_   = globalproblem_->ImmersedMethodParams().get<std::string>("APPLY_VEL_RELAX")   == "globally";
    forcerelax_         = globalproblem_->ImmersedMethodParams().get<double>("FORCE_RELAX");
    velrelax_           = globalproblem_->ImmersedMethodParams().get<double>("VEL_RELAX");

    if(globalproblem_->CellMigrationParams().sublist("FLOW INTERACTION MODULE").get<std::string>("COUPALGO") == "iter_stagg_fixed_rel_param")
    {
      coupalgo_= INPAR::IMMERSED::fixed;
      if(myrank_==0)
        std::cout<<"\n"
                 " Relax Force globally = "<<relaxforceglobally_<<" with relaxation parameter = "<<forcerelax_<<"\n"
                 " Relax Vel   globally = "<<relaxvelglobally_  <<" with relaxation parameter = "<<velrelax_<<"\n"<<std::endl;
    }
    else if(globalproblem_->CellMigrationParams().sublist("FLOW INTERACTION MODULE").get<std::string>("COUPALGO") == "iter_stagg_AITKEN_rel_param")
    {
      coupalgo_ = INPAR::IMMERSED::aitken;
      if(myrank_==0)
        std::cout<<"\n Using AITKEN relaxation parameter. "<<std::endl;
    }
    else
      dserror("Unknown definition of COUPALGO in FSI DYNAMIC section for Immersed FSI.");


  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::ReadRestart(int step)
{
  cellstructure_->ReadRestart(step);
  poroscatra_subproblem_->ReadRestart(step);
  SetTimeStep(poroscatra_subproblem_->PoroField()->Time(),step);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::CalcFluidTractionOnStructure()
{
  Teuchos::ParameterList params;
  params.set<std::string>("action","calc_fluid_traction");
  params.set<std::string>("backgrddisname","porofluid");
  params.set<std::string>("immerseddisname","cell");

  backgroundfluiddis_->SetState(0,"velnp",poroscatra_subproblem_->FluidField()->Velnp());
  backgroundfluiddis_->SetState(0,"veln", poroscatra_subproblem_->FluidField()->Veln());

  double normofstructbdrytraction=-1234.0;

  DRT::AssembleStrategy cell_fld_bdry_strategy(
      0,              // struct dofset for row
      0,              // struct dofset for column
      Teuchos::null,  // matrix 1
      Teuchos::null,  //
      cell_bdry_traction_, // vector 1
      Teuchos::null,  //
      Teuchos::null   //
  );

  if(myrank_ == 0)
  {
    std::cout<<"################################################################################################"<<std::endl;
    std::cout<<"###   Interpolate fluid stresses to structural surface and calculate tractions                  "<<std::endl;

  }
  EvaluateInterpolationCondition( immerseddis_, params, cell_fld_bdry_strategy, "IMMERSEDCoupling", -1 );

  cell_bdry_traction_->Norm2(&normofstructbdrytraction);
  if(myrank_ == 0)
  {
    std::cout<<"###   Norm of Boundary Fluid Traction:   "<<std::setprecision(10)<<normofstructbdrytraction<<std::endl;
    std::cout<<"################################################################################################"<<std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedFlowCellInteraction::ResetImmersedInformation()
{
  Teuchos::ParameterList params;
  params.set<int>("action",FLD::reset_immersed_ele);
  params.set<int>("Physical Type", poroscatra_subproblem_->FluidField()->PhysicalType());
  params.set<int>("intpoints_fluid_bound",degree_gp_fluid_bound_);
  backgroundfluiddis_->Evaluate(params);

  return;
}
