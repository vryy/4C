/*!----------------------------------------------------------------------
\file immersed_partitioned_cellmigration.cpp

\brief partitioned immersed cell migration algorithm

<pre>
Maintainers: Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "immersed_partitioned_cellmigration.H"

#include "../linalg/linalg_utils.H"

#include "../drt_poroelast/poro_scatra_base.H"
#include "../drt_poroelast/poroelast_utils_setup.H"

#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"
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


IMMERSED::ImmersedPartitionedCellMigration::ImmersedPartitionedCellMigration(const Epetra_Comm& comm)
  : ImmersedPartitioned(comm)
{
  // important variables for parallel simulations
  myrank_  = comm.MyPID();
  numproc_ = comm.NumProc();

  // get pointer to global problem
  globalproblem_ = DRT::Problem::Instance();

  // get pointer to discretizations
  backgroundfluiddis_     = globalproblem_->GetDis("porofluid");
  backgroundstructuredis_ = globalproblem_->GetDis("structure");
  immerseddis_            = globalproblem_->GetDis("cell");
  scatradis_              = globalproblem_->GetDis("scatra");

  // counter for continued (unconverged) steps
  continued_steps_=0;

  // decide whether multiple structural bodies or not
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
    dserror("set MIGRATIONTYPE to 'ameboid' or 'proteolytic' in --CELL DYNAMICS sectionin your .dat file.");

  // initialize segregation variables
  segregationconstant_=globalproblem_->CellMigrationParams().get<double>("SEGREGATION_CONST");
  segregationtype_=DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(),"SEGREGATION");
  segregationby_=DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(),"SEGREGATION_BY");

  // setup the relaxation parameters
  SetupRelaxation();

  // build field  cell structure
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> cellstructure =
        Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(globalproblem_->CellMigrationParams(), const_cast<Teuchos::ParameterList&>(globalproblem_->StructuralDynamicParams()), immerseddis_));
    cellstructure_ = Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapperImmersed>(cellstructure->StructureField());
    if(cellstructure_==Teuchos::null)
      dserror("dynamic cast from Structure to FSIStructureWrapperImmersed failed");
    if(myrank_==0)
      std::cout<<"Created Field Cell Structure ..."<<std::endl;

  // create instance of poroelast subproblem
  //poroelast_subproblem_ = Teuchos::rcp(new POROELAST::Monolithic(comm, globalproblem_->PoroelastDynamicParams()));
  poroelast_subproblem_ = POROELAST::UTILS::CreatePoroScatraAlgorithm(globalproblem_->PoroScatraControlParams(),Comm());
  // set pointer to poro fpsi structure
  porostructure_ = poroelast_subproblem_->PoroField()->StructureField();
  // setup of poro monolithic algorithm
  SetupPoroAlgorithm();
  if(myrank_==0)
    std::cout<<" Created Field Poroelast ... \n"<<std::endl;

  // vector of fluid stresses interpolated to cell bdry int points integrated over cell surface
  cell_bdry_traction_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()),true));
  // vector of penalty traction on cell bdry int points integrated over cell surface
  cell_penalty_traction_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()),true));
  // vector of DELTA of penalty traction on cell bdry int points integrated over cell surface
  cell_delta_penalty_traction_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()),true));
  // vector of DELTA of penalty traction on cell bdry int points integrated over cell surface of previous iteration
  cell_delta_penalty_traction_old_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()),true));
  // gap integrated over surface
  penalty_gap_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()),true));
  // current nodal normals
  curr_nodal_normals_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()),true));
  // vector with fluid velocities interpolated from structure
  porofluid_artificial_velocity_ = Teuchos::rcp(new Epetra_Vector(*(poroelast_subproblem_->FluidField()->DofRowMap()),true));
  // vector with fluid velocities interpolated from structure
  porofluid_artificial_source_ = Teuchos::rcp(new Epetra_Vector(*(poroelast_subproblem_->FluidField()->DofRowMap()),true));
  // vector with fluid velocities interpolated from structure
  poroscatra_segregated_phi_ = Teuchos::rcp(new Epetra_Vector(*(poroelast_subproblem_->ScaTraField()->DofRowMap(0)),true));

  // construct 3D search tree for fluid domain
  // initialized in SetupSBackgroundDiscretization()
  fluid_SearchTree_ = Teuchos::rcp(new GEO::SearchTree(5));

  // construct 3D search tree for cell domain
  // initialized in SetupImmersedDiscretization()
  cell_SearchTree_ = Teuchos::rcp(new GEO::SearchTree(5));

  // get coupling variable
  displacementcoupling_ = globalproblem_->ImmersedMethodParams().sublist("PARTITIONED SOLVER").get<std::string>("COUPVARIABLE") == "Displacement";
  if(displacementcoupling_ and myrank_==0)
    std::cout<<" Coupling variable for partitioned scheme :  Displacements "<<std::endl;
  else if (!displacementcoupling_ and myrank_==0)
    std::cout<<" Coupling variable for partitioned scheme :  Force "<<std::endl;

  // penalty constraint enforcement parameters
  penalty_start_ = globalproblem_->CellMigrationParams().get<double>("PENALTY_START");
  penalty_init_  = globalproblem_->CellMigrationParams().get<double>("PENALTY_INIT");

  // construct immersed exchange manager. singleton class that makes immersed variables comfortably accessible from everywhere in the code
  exchange_manager_ = DRT::ImmersedFieldExchangeManager::Instance();
  exchange_manager_->InitializePorosityAtGPMap();
  exchange_manager_->SetPointerToECMPenaltyTraction(cell_penalty_traction_);
  exchange_manager_->SetPointerToGap(penalty_gap_);
  exchange_manager_->SetPointerToCurrentNodalNormals(curr_nodal_normals_);
  exchange_manager_->SetIsInitialized(true);

  // initialize the parameter initielize_cell_ which determines whether or not first time step is pre-simulation
  initialize_cell_=DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(),"INITIALIZE_CELL");
  timestep_ = Dt();
  if(initialize_cell_)
    initial_timestep_=globalproblem_->CellMigrationParams().get<double>("INITIAL_TIMESTEP");
  else
    initial_timestep_=timestep_;

  initialization_steps_=globalproblem_->CellMigrationParams().get<int>("INITIALIZATION_STEPS");

  ecm_interaction_=DRT::INPUT::IntegralValue<int>(globalproblem_->CellMigrationParams(),"ECM_INTERACTION");

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
    dserror("Invalid value for parameter NUM_GP_FLUID_BOUND (valid parameters are 8, 64, 125, 343, 729 and 1000).");

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::CouplingOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  // reinitialize the transfer vectors
  porofluid_artificial_velocity_->PutScalar(0.0);
  porofluid_artificial_source_->PutScalar(0.0);
  cell_bdry_traction_->PutScalar(0.0);
  cell_delta_penalty_traction_->PutScalar(0.0);
  poroscatra_segregated_phi_->PutScalar(0.0);
  penalty_gap_->PutScalar(0.0);

  // reset element and node information about immersed method
  Teuchos::ParameterList params;
  params.set<int>("action",FLD::reset_immersed_ele);
  params.set<int>("Physical Type", poroelast_subproblem_->FluidField()->PhysicalType());
  params.set<int>("intpoints_fluid_bound",degree_gp_fluid_bound_);
  backgroundfluiddis_->Evaluate(params);


  if (displacementcoupling_) // DISPLACEMENT COUPLING
  {
    const Teuchos::RCP<Epetra_Vector> idispn = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL BackgroundOp
    ////////////////////
    PrepareBackgroundOp();
    if(!initialize_cell_ or Step()>initialization_steps_)
    {
      BackgroundOp(porofluid_artificial_velocity_, fillFlag);

    // reset element and node information about immersed method
    params.set<int>("action",FLD::reset_immersed_ele);
    poroelast_subproblem_->FluidField()->Discretization()->Evaluate(params);

    }

    ////////////////////
    // CALL ImmersedOp
    ////////////////////
    PrepareImmersedOp();
    Teuchos::RCP<Epetra_Vector> idispnp =
    ImmersedOp(cell_bdry_traction_, fillFlag);

    int err = F.Update(1.0, *idispnp, -1.0, *idispn, 0.0);
    if(err != 0)
      dserror("Vector update of Coupling-residual returned err=%d",err);

  }
  else if(!displacementcoupling_) // FORCE COUPLING
  {

    // get the initial guess
    const Teuchos::RCP<Epetra_Vector> iforcen = Teuchos::rcp(new Epetra_Vector(x));


    ////////////////////
    // CALL ImmersedOp
    ////////////////////
    PrepareImmersedOp();
    ImmersedOp(cell_bdry_traction_, fillFlag);

    ////////////////////
    // CALL BackgroundOp
    ////////////////////
    PrepareBackgroundOp(); // determine elements cut by the boundary
    if(!initialize_cell_ or Step()>initialization_steps_)
    {
      BackgroundOp(porofluid_artificial_velocity_, fillFlag);

      // reset element and node information about immersed method
      params.set<int>("action",FLD::reset_immersed_ele);
      params.set<int>("intpoints_fluid_bound",degree_gp_fluid_bound_);
      backgroundfluiddis_->Evaluate(params);
    }

    ///////////////////////////////////////////////////
    // CALL ImmersedOp again to recalculate tractions (only tractions due to ImmersedOp( , User))
    //////////////////////////////////////////////////
    const Teuchos::RCP<Epetra_Vector> iforcenp = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()),true));
    ImmersedOp(iforcenp, User);

    int err = F.Update(1.0, *(cellstructure_->Interface()->ExtractIMMERSEDCondVector(iforcenp)), -1.0, *iforcen, 0.0);

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
void IMMERSED::ImmersedPartitionedCellMigration::BackgroundOp(Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values,
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


    // relaxation
    if(displacementcoupling_ and abs(velrelax_-1.0)>1e-13 and relaxvelglobally_ and coupalgo_ == INPAR::IMMERSED::fixed)
    {
      std::cout<<"Apply relaxation parameter omega = "<<velrelax_<<" on full artificial velocity field."<<std::endl;
      backgrd_dirichlet_values->Scale(velrelax_);
    }
    else if(displacementcoupling_ and abs(velrelax_-1.0)>1e-13 and relaxvelglobally_ and coupalgo_ == INPAR::IMMERSED::aitken)
    {
      double omega_aitken = Teuchos::rcp_dynamic_cast<NOX::FSI::AitkenFactory>(aitkenfactory_)->GetAitken()->getOmega();
      std::cout<<"Apply AITKEN relaxation parameter omega = "<<omega_aitken<<" on full artificial velocity field."<<std::endl;

      backgrd_dirichlet_values->Scale(omega_aitken);
    }
    else if(displacementcoupling_ and relaxvelglobally_ == false)
      std::cout<<"Selective velocity relaxation not implemented yet. Parameter VEL_RELAX = "<<velrelax_<<" has no effect."<<std::endl;


    if(myrank_ == 0)
    {
      std::cout<<"################################################################################################"<<std::endl;
      std::cout<<"###   Interpolate transfer quantities from immersed elements which overlap the "<<backgroundfluiddis_->Name()<<" nodes ..."<<std::endl;
    }

    if(curr_subset_of_backgrounddis_.empty()==false)
      EvaluateWithInternalCommunication(params,
                                        backgroundfluiddis_,
                                        &fluid_vol_strategy,
                                        curr_subset_of_backgrounddis_,
                                        cell_SearchTree_, currpositions_cell_,
                                        (int)FLD::interpolate_velocity_to_given_point,
                                        false);

    /////////////////////////////////////////////////////////////////////////////////////
    //  Apply Dirichlet values to background fluid
    //
    // dofsetnum of backgroundfluid is 0 for backgroundfluiddis in poroelast_subproblem_
    BuildImmersedDirichMap(backgroundfluiddis_, dbcmap_immersed_, poroelast_subproblem_->FluidField()->GetDBCMapExtractor()->CondMap(),0);
    poroelast_subproblem_->FluidField()->AddDirichCond(dbcmap_immersed_);
    // apply immersed dirichlets to porofluid field
    DoImmersedDirichletCond(poroelast_subproblem_->FluidField()->WriteAccessVelnp(),backgrd_dirichlet_values, dbcmap_immersed_);

//    /////////////////////////////////////////////////////////////////////////////////////
//    //  Apply Source values to rhs of background fluid
//    //
//    backgroundstructuredis_->SetState(0,"dispnp",porostructure_->Dispnp());
//    backgroundfluiddis_->SetState(0,"velnp",poroelast_subproblem_->FluidField()->Velnp());
//
//    if(curr_subset_of_backgrounddis_.empty()==false)
//      EvaluateWithInternalCommunication(backgroundfluiddis_,
//                                        &fluid_vol_strategy,
//                                        curr_subset_of_backgrounddis_,
//                                        cell_SearchTree_, currpositions_cell_,
//                                        (int)FLD::calc_artificial_velocity_divergence,
//                                        false);
//
//    ExtractPressurePart();
//    double normofsource=-1234.0;
//    porofluid_artificial_source_->Norm2(&normofsource);
//    poroelast_subproblem_->FluidField()->AddContributionToExternalLoads(porofluid_artificial_source_);

    /////////////////////////////////////////////////////////////////////////////////////
    //  Apply Dirichlet values to background scatra field
    //
    if(migrationtype_==INPAR::CELL::cell_migration_proteolytic and segregationby_==INPAR::CELL::segregation_by_dirichlet and segregationtype_==INPAR::CELL::segregation_volumetric)
    {
      // only constant segregation so far
      poroscatra_segregated_phi_->PutScalar(segregationconstant_);
      BuildImmersedScaTraDirichMap(backgroundfluiddis_, scatradis_, dbcmap_immersed_scatra_, poroelast_subproblem_->ScaTraField()->DirichMaps()->CondMap(),0);
      poroelast_subproblem_->ScaTraField()->AddDirichCond(dbcmap_immersed_scatra_);
      DoImmersedDirichletCond(poroelast_subproblem_->ScaTraField()->Phinp(),poroscatra_segregated_phi_, dbcmap_immersed_scatra_);
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

      // set poro solution in scatra
      poroelast_subproblem_->SetPoroSolution();

      Teuchos::ParameterList sparams_struct;

      scatradis_->AddMultiVectorToParameterList(sparams_struct,
                                                "dispnp",
                                                poroelast_subproblem_->ScaTraField()->Disp()
                                               );

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
                                                curr_subset_of_backgrounddis_,
                                                cell_SearchTree_,
                                                currpositions_cell_,
                                                sparams_struct,
                                                evaluateonlyboundary);

      poroelast_subproblem_->ScaTraField()->AddContributionToRHS(poroscatra_segregated_phi_);

    }
    else
    {
      if(migrationtype_==INPAR::CELL::cell_migration_proteolytic)
        dserror("combination of SEGREGATION parameters is not implemented. Fix your input file.");
    }


    // rebuild the combined dbcmap
    poroelast_subproblem_->BuildCombinedDBCMap();

    double normofvelocities = -1234.0;
    poroelast_subproblem_->FluidField()->ExtractVelocityPart(backgrd_dirichlet_values)->Norm2(&normofvelocities);

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
    exchange_manager_->SetCheckCounter(0);
    //exchange_manager_->SetNumIsImmersedBoundary(0);
    poroelast_subproblem_->Solve();
    //std::cout<<"Matched "<<exchange_manager_->GetCheckCounter()<<" integration points in "<<exchange_manager_->GetNumIsImmersedBoundary()<<" IsImmersedBoundary() elements on PROC "<<myrank_<<". "<<std::endl;

    // remove immersed dirichlets from dbcmap of fluid (and scatra) (may be different in next iteration)
    poroelast_subproblem_->FluidField() ->RemoveDirichCond(dbcmap_immersed_);
    if(migrationtype_==INPAR::CELL::cell_migration_proteolytic and segregationby_==INPAR::CELL::segregation_by_dirichlet)
      poroelast_subproblem_->ScaTraField()->RemoveDirichCond(dbcmap_immersed_scatra_);

    // rebuild the combined dbcmap
    poroelast_subproblem_->BuildCombinedDBCMap();

  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedCellMigration::ImmersedOp(Teuchos::RCP<Epetra_Vector> bdry_traction,
                                                       const FillType fillFlag)
{
  IMMERSED::ImmersedPartitioned::ImmersedOp(bdry_traction,fillFlag);

  if(fillFlag==User)
  {
    Teuchos::ParameterList params;
    params.set<std::string>("action","calc_fluid_traction");
    params.set<std::string>("backgrddisname","porofluid");
    params.set<std::string>("immerseddisname","cell");

    backgroundfluiddis_->SetState(0,"velnp",poroelast_subproblem_->FluidField()->Velnp());
    backgroundfluiddis_->SetState(0,"veln", poroelast_subproblem_->FluidField()->Veln());

    DRT::AssembleStrategy struct_bdry_strategy(
        0,              // struct dofset for row
        0,              // struct dofset for column
        Teuchos::null,  // matrix 1
        Teuchos::null,  //
        bdry_traction , // vector 1
        Teuchos::null,  //
        Teuchos::null   //
    );
    if(myrank_ == 0)
    {
      std::cout<<"################################################################################################"<<std::endl;
      std::cout<<"###   Calculate Boundary Tractions on Cell for Initial Guess or Convergence Check ...      "<<std::endl;

    }
    EvaluateInterpolationCondition( immerseddis_, params, struct_bdry_strategy,"IMMERSEDCoupling", -1 );
    bdry_traction->Update(1.0,*cell_penalty_traction_,1.0);

    double normorstructbdrytraction;
    bdry_traction->Norm2(&normorstructbdrytraction);
    if(myrank_ == 0 and ecm_interaction_)
    {
      std::cout<<"###   Norm of Boundary Fluid+ECM Traction:   "<<std::setprecision(10)<<normorstructbdrytraction<<std::endl;
      std::cout<<"################################################################################################"<<std::endl;
    }
    else if(myrank_ == 0 and !ecm_interaction_)
    {
      std::cout<<"###   Norm of Boundary Fluid Traction:   "<<std::setprecision(10)<<normorstructbdrytraction<<std::endl;
      std::cout<<"################################################################################################"<<std::endl;
    }

    return Teuchos::null;
  }
  else
  {
      Teuchos::ParameterList params;
      params.set<std::string>("action","calc_fluid_traction");
      params.set<std::string>("backgrddisname","porofluid");
      params.set<std::string>("immerseddisname","cell");

      backgroundfluiddis_->SetState(0,"velnp",poroelast_subproblem_->FluidField()->Velnp());
      backgroundfluiddis_->SetState(0,"veln", poroelast_subproblem_->FluidField()->Veln());

      double normofstructbdrytraction=-1234.0;
      double normofdeltastructbdrytraction=-1234.0;


      if(!initialize_cell_) // not during cell initialization
      {

        DRT::AssembleStrategy cell_fld_bdry_strategy(
            0,              // struct dofset for row
            0,              // struct dofset for column
            Teuchos::null,  // matrix 1
            Teuchos::null,  //
            bdry_traction, // vector 1
            Teuchos::null,  //
            Teuchos::null   //
        );

        if(myrank_ == 0)
        {
          std::cout<<"################################################################################################"<<std::endl;
          std::cout<<"###   Interpolate fluid stresses to structural surface and calculate tractions                  "<<std::endl;

        }
        EvaluateInterpolationCondition( immerseddis_, params, cell_fld_bdry_strategy, "IMMERSEDCoupling", -1 );

        bdry_traction->Norm2(&normofstructbdrytraction);
        if(myrank_ == 0)
        {
          std::cout<<"###   Norm of Boundary Fluid Traction:   "<<std::setprecision(10)<<normofstructbdrytraction<<std::endl;
          std::cout<<"################################################################################################"<<std::endl;
        }
      } // if not during cell initialization


    ///////////////////////////////////////////7
    //  Cell - ECM Interaction Part
    ////////////////////////////////////////////
    if(ecm_interaction_)
    {
//      DRT::AssembleStrategy cell_str_nodal_normal_strategy(
//          0,              // struct dofset for row
//          0,              // struct dofset for column
//          Teuchos::null,  // matrix 1
//          Teuchos::null,  //
//          curr_nodal_normals_, // vector 1
//          Teuchos::null,  // vector 2
//          Teuchos::null   //
//      );
//
//      params.set<std::string>("action","calc_cur_nodal_normals");
//      EvaluateInterpolationCondition( immerseddis_, params, cell_str_nodal_normal_strategy, "IMMERSEDCoupling", -1 );

      DRT::AssembleStrategy cell_str_bdry_strategy(
          0,              // struct dofset for row
          0,              // struct dofset for column
          Teuchos::null,  // matrix 1
          Teuchos::null,  //
          cell_delta_penalty_traction_ , // vector 1
          penalty_gap_,                  // vector 2
          Teuchos::null   //
      );
      if(myrank_ == 0)
      {
        std::cout<<"################################################################################################"<<std::endl;
        std::cout<<"###   Calculate ECM stresses on structural surface and calculate tractions                  "<<std::endl;

      }

      exchange_manager_->SetGapMax(-1234.0);
      exchange_manager_->SetVoidMax(-1234.0);
      exchange_manager_->SetGapMin(1234.0);
      exchange_manager_->SetVoidMin(1234.0);
      exchange_manager_->SetDeltaPorosityMax(0.0);

      double penalty=penalty_start_;
      params.set<std::string>("action","calc_ecm_traction");
      params.set<std::string>("backgrddisname","structure");

      if(!initialize_cell_)
        penalty=penalty_start_;//*sqrt(((IterationCounter())[0]));
      else
        penalty=penalty_init_;

      params.set<double>("penalty",penalty);
      params.set<int>("during_init",initialize_cell_);
      EvaluateInterpolationCondition( immerseddis_, params, cell_str_bdry_strategy, "IMMERSEDCoupling", -1 );

//      if(!initialize_cell_ and (IterationCounter()[0])>1)
//      {
//        // AITKEN relaxation for ecm interaction force
//        double top = 1.0;
//        double den = 1.0;
//        Teuchos::RCP<Epetra_Vector> temp = Teuchos::rcp(new Epetra_Vector(*cell_delta_penalty_traction_old_));
//        int err = temp->Update(-1.0,*cell_delta_penalty_traction_,1.0);
//        if(err != 0)
//          dserror("Epetra Update returned error code err=%d",err);
//
//          temp->Norm2(&den);
//
//        err = cell_delta_penalty_traction_old_->Dot(*temp,&top);
//        if(err != 0)
//          dserror("Epetra Dot Product returned error code err=%d",err);
//
//        if(abs(den)<1.0e-12)
//          dserror("denominator too small den=%f",den);
//        nu_ = nu_ + (nu_ - 1.0)*top/(den*den);
//
//        cell_delta_penalty_traction_old_->Update(1.0,*cell_delta_penalty_traction_,0.0);
//      }

      if(!initialize_cell_)
      {
          //DeleteForceAtLeadingEdge(cell_delta_penalty_traction_);
          //DeleteForceAtLeadingEdge(penalty_gap_);
      }

      // Apply relaxed ecm force increment
      //PreventInterfacePullingTraction(cell_delta_penalty_traction_,cell_penalty_traction_);
      cell_penalty_traction_->Update((1.0-nu_),*cell_delta_penalty_traction_,0.0); // cell_delta_penalty_traction contains now full force, thus this*0.0
      cell_delta_penalty_traction_->Norm2(&normofdeltastructbdrytraction);

      bdry_traction->Update(1.0,*cell_penalty_traction_,1.0);
      normofstructbdrytraction = 0.0;
      cell_penalty_traction_->Norm2(&normofstructbdrytraction);
      if(myrank_ == 0)
      {
        std::cout<<"###   Penalty Parameter:  "<<penalty<<std::endl;
        //std::cout<<"###   AITKEN Step Parameter:  "<<std::setprecision(6)<<1.0-nu_<<std::endl;
        std::cout<<"###   Max delta_phi:  "<<std::setprecision(6)<<exchange_manager_->GetDeltaPorosityMax()<<std::endl;
        std::cout<<"###   MAX gap  = "<<std::setprecision(6)<<exchange_manager_->GetGapMax()<<"   MIN gap  = "<<std::setprecision(6)<<exchange_manager_->GetGapMin()<<std::endl;
        std::cout<<"###   Max gap origin:  "<<std::setprecision(6)<<*(exchange_manager_->GetMaxGapSpacePoint())<<std::endl;
        std::cout<<"###   Max gap end point:  "<<*(exchange_manager_->GetMaxGapSearchDirection())<<std::endl;
        std::cout<<"###   Min gap origin:  "<<std::setprecision(6)<<*(exchange_manager_->GetMinGapSpacePoint())<<std::endl;
        std::cout<<"###   Min gap end point:  "<<*(exchange_manager_->GetMinGapSearchDirection())<<std::endl;
        std::cout<<"###   MAX void = "<<std::setprecision(6)<<exchange_manager_->GetVoidMax()<<"   MIN void = "<<std::setprecision(6)<<exchange_manager_->GetVoidMin()<<std::endl;
        std::cout<<"###   Norm of Boundary ECM interaction force vector:    "<<std::setprecision(13)<<normofstructbdrytraction<<std::endl;
        std::cout<<"###   Norm of Boundary ECM interaction force increment: "<<std::setprecision(13)<<normofdeltastructbdrytraction<<std::endl;
        std::cout<<"################################################################################################"<<std::endl;
      }
    } // if cell-ecm interaction is turned on


    // relaxation
    if(!displacementcoupling_ and abs(velrelax_-1.0)>1e-13 and relaxvelglobally_ and coupalgo_ == INPAR::IMMERSED::fixed)
    {
      std::cout<<"Apply relaxation parameter omega = "<<velrelax_<<" on whole structural force interface."<<std::endl;
      bdry_traction->Scale(velrelax_);
    }
    else if(!displacementcoupling_ and abs(velrelax_-1.0)>1e-13 and relaxvelglobally_ and coupalgo_ == INPAR::IMMERSED::aitken)
    {
      double omega_aitken = Teuchos::rcp_dynamic_cast<NOX::FSI::AitkenFactory>(aitkenfactory_)->GetAitken()->getOmega();
      std::cout<<"Apply AITKEN relaxation parameter omega = "<<omega_aitken<<" on whole structural force interface."<<std::endl;

      bdry_traction->Scale(omega_aitken);
    }
    else if(!displacementcoupling_ and relaxvelglobally_ == false)
      std::cout<<"Selective velocity relaxation not implemented yet. Parameter VEL_RELAX = "<<velrelax_<<" has no effect."<<std::endl;

    // prescribe neumann values at structural boundary dofs
    cellstructure_->ApplyImmersedInterfaceForces(Teuchos::null,cellstructure_->Interface()->ExtractIMMERSEDCondVector(*bdry_traction));

    // solve
    cellstructure_->Solve();

    return cellstructure_->ExtractImmersedInterfaceDispnp();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  if(myrank_==0)
    std::cout<<"Cell Predictor: "<<std::endl;
  cellstructure_->PrepareTimeStep();
  if(myrank_==0)
    std::cout<<"Poro Predictor: "<<std::endl;
  poroelast_subproblem_->PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedCellMigration::InitialGuess()
{
  if(displacementcoupling_)
    return cellstructure_->PredictImmersedInterfaceDispnp();
  else // FORCE COUPLING
  {
    immerseddis_->SetState(0,"displacement",cellstructure_->Dispnp());
    cell_bdry_traction_old_ = Teuchos::rcp(new Epetra_Vector(*(cellstructure_->DofRowMap()),true));
    ImmersedOp(cell_bdry_traction_old_, User);

    return cellstructure_->Interface()->ExtractIMMERSEDCondVector(cell_bdry_traction_old_);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::BuildImmersedDirichMap(Teuchos::RCP<DRT::Discretization> dis,
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

          // include also pressure dof if node does not belong to a boundary background element
          //if((nodes[inode]->IsBoundaryImmersed())==0)
          //  mydirichdofs.push_back(dofs[3]);
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
void IMMERSED::ImmersedPartitionedCellMigration::BuildImmersedScaTraDirichMap(Teuchos::RCP<DRT::Discretization> immersedinfodis,
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
void IMMERSED::ImmersedPartitionedCellMigration::DoImmersedDirichletCond(Teuchos::RCP<Epetra_Vector> statevector, Teuchos::RCP<Epetra_Vector> dirichvals, Teuchos::RCP<Epetra_Map> dbcmap)
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
void IMMERSED::ImmersedPartitionedCellMigration::SetupImmersedDiscretization()
{
  // ghost structure on each proc (for search algorithm)
  if(numproc_ > 1)
  {
    // fill complete inside
    CreateGhosting(immerseddis_);
  }
  else
  {
    // fill complete to incorporate changes due to ghosting and build geometries
    immerseddis_->FillComplete();
  }

  // find positions of the immersed discretization on proc
  std::map<int,LINALG::Matrix<3,1> > my_currpositions_cell;
  for (int lid = 0; lid < immerseddis_->NumMyRowNodes(); ++lid)
  {
    const DRT::Node* node = immerseddis_->lRowNode(lid);
    LINALG::Matrix<3,1> currpos;

    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];

    my_currpositions_cell[node->Id()] = currpos;
  }
  // Communicate local currpositions:
  // map with current structural positions should be same on all procs
  // to make use of the advantages of ghosting the structure redundantly
  // on all procs.
  int procs[numproc_];
  for(int i=0;i<numproc_;i++)
    procs[i]=i;
  LINALG::Gather<int,LINALG::Matrix<3,1> >(my_currpositions_cell,currpositions_cell_,numproc_,&procs[0],Comm());

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3,2> rootBox2 = GEO::getXAABBofDis(*immerseddis_,currpositions_cell_);
  cell_SearchTree_->initializeTree(rootBox2,*immerseddis_,GEO::TreeType(GEO::OCTTREE));

  if(myrank_==0)
    std::cout<<"\n Build Cell SearchTree ... "<<std::endl;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::SetupBackgroundDiscretization()
{
  // find positions of the background fluid discretization
  for (int lid = 0; lid < backgroundfluiddis_->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = backgroundfluiddis_->lColNode(lid);
    LINALG::Matrix<3,1> currpos;

    currpos(0) = node->X()[0];
    currpos(1) = node->X()[1];
    currpos(2) = node->X()[2];

    currpositions_ECM_[node->Id()] = currpos;
  }

  // find the bounding box of the elements and initialize the search tree
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*backgroundfluiddis_,currpositions_ECM_);
  fluid_SearchTree_->initializeTree(rootBox,*backgroundfluiddis_,GEO::TreeType(GEO::OCTTREE));

  if(myrank_==0)
    std::cout<<"\n Build Fluid SearchTree ... "<<std::endl;
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::PrepareImmersedOp()
{
  backgroundstructuredis_->SetState(0,"displacement",porostructure_->Dispnp());
  immerseddis_->SetState(0,"velnp",cellstructure_->Velnp());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::PrepareBackgroundOp()
{
  TEUCHOS_FUNC_TIME_MONITOR("IMMERSED::PrepareBackgroundOp()");

  immerseddis_->SetState(0,"displacement",cellstructure_->Dispnp());
  immerseddis_->SetState(0,"displacement_old",cellstructure_->Dispn());
  immerseddis_->SetState(0,"velocity_old",cellstructure_->Veln());
  immerseddis_->SetState(0,"velocity",cellstructure_->Velnp());
  backgroundfluiddis_->SetState(0,"veln",poroelast_subproblem_->FluidField()->Veln());

  double structsearchradiusfac = DRT::Problem::Instance()->ImmersedMethodParams().get<double>("STRCT_SRCHRADIUS_FAC");

  //
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
  // Communicate local currpositions:
  // map with current structural positions should be same on all procs
  // to make use of the advantages of ghosting the structure redundantly
  // on all procs.
  int procs[Comm().NumProc()];
  for(int i=0;i<Comm().NumProc();i++)
    procs[i]=i;
  LINALG::Gather<int,LINALG::Matrix<3,1> >(my_currpositions_cell,currpositions_cell_,Comm().NumProc(),&procs[0],Comm());

//  DEBUG output
//  std::cout<<"Proc "<<Comm().MyPID()<<": my_curr_pos.size()="<<my_currpositions_cell.size()<<std::endl;
//  std::cout<<"Proc "<<Comm().MyPID()<<":    curr_pos.size()="<<currpositions_struct_.size()<<std::endl;

  if (multicellmigration_ == false)
  {
    // get bounding box of current configuration of structural dis
    const LINALG::Matrix<3,2> structBox = GEO::getXAABBofDis(*immerseddis_,currpositions_cell_);
    double max_radius = sqrt(pow(structBox(0,0)-structBox(0,1),2)+pow(structBox(1,0)-structBox(1,1),2)+pow(structBox(2,0)-structBox(2,1),2));
    // search for background elements within a certain radius around the center of the immersed bounding box
    LINALG::Matrix<3,1> boundingboxcenter;
    boundingboxcenter(0) = structBox(0,0)+(structBox(0,1)-structBox(0,0))*0.5;
    boundingboxcenter(1) = structBox(1,0)+(structBox(1,1)-structBox(1,0))*0.5;
    boundingboxcenter(2) = structBox(2,0)+(structBox(2,1)-structBox(2,0))*0.5;

# ifdef DEBUG
    std::cout<<"Bounding Box of Structure: "<<structBox<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Bounding Box Center of Structure: "<<boundingboxcenter<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Search Radius Around Center: "<<structsearchradiusfac*max_radius<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Length of Dispnp()="<<celldisplacements->MyLength()<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Size of currpositions_struct_="<<currpositions_cell_.size()<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"My DofRowMap Size="<<immerseddis_->DofRowMap()->NumMyElements()<<"  My DofColMap Size="<<immerseddis_->DofColMap()->NumMyElements()<<" on PROC "<<Comm().MyPID()<<std::endl;
    std::cout<<"Dis Structure NumColEles: "<<immerseddis_->NumMyColElements()<<" on PROC "<<Comm().MyPID()<<std::endl;
# endif

    curr_subset_of_backgrounddis_ = fluid_SearchTree_->searchElementsInRadius(*backgroundfluiddis_,currpositions_ECM_,boundingboxcenter,structsearchradiusfac*max_radius,0);

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
void IMMERSED::ImmersedPartitionedCellMigration::SetupPoroAlgorithm()
{
  poroelast_subproblem_->SetupSystem();
  //poroelast_subproblem_->SetupSolver(); // done in CreatePoroScatraAlgorithm
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::PrepareOutput()
{
  poroelast_subproblem_->PrepareOutput();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::Output()
{
  cellstructure_->Output();
  poroelast_subproblem_->Output();

  if(myrank_==0)
    std::cout<<" Number of unconverged steps: "<<continued_steps_<<"\n"<<std::endl;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::Update()
{
  cellstructure_->Update();
  poroelast_subproblem_->Update();
  exchange_manager_->UpdatePorosityAtGPOldTimestep();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::SetupRelaxation()
{
  // initialize some relaxation related member variables
    relaxforceglobally_ = globalproblem_->ImmersedMethodParams().get<std::string>("APPLY_FORCE_RELAX") == "globally";
    relaxvelglobally_   = globalproblem_->ImmersedMethodParams().get<std::string>("APPLY_VEL_RELAX")   == "globally";
    forcerelax_         = globalproblem_->ImmersedMethodParams().get<double>("FORCE_RELAX");
    velrelax_           = globalproblem_->ImmersedMethodParams().get<double>("VEL_RELAX");

    if(globalproblem_->CellMigrationParams().get<std::string>("COUPALGO") == "iter_stagg_fixed_rel_param")
    {
      coupalgo_= INPAR::IMMERSED::fixed;
      if(myrank_==0)
        std::cout<<"\n"
                 " Relax Force globally = "<<relaxforceglobally_<<" with relaxation parameter = "<<forcerelax_<<"\n"
                 " Relax Vel   globally = "<<relaxvelglobally_  <<" with relaxation parameter = "<<velrelax_<<"\n"<<std::endl;
    }
    else if(globalproblem_->CellMigrationParams().get<std::string>("COUPALGO") == "iter_stagg_AITKEN_rel_param")
    {
      coupalgo_ = INPAR::IMMERSED::aitken;
      if(myrank_==0)
        std::cout<<"\n Using AITKEN relaxation parameter. "<<std::endl;
    }
    else
      dserror("Unknown definition of COUPALGO in FSI DYNAMIC section for Immersed FSI.");

    // internal AITKEN parameter used for relaxation of ecm force increment
    nu_=0.0;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::ExtractPressurePart()
{
  double* values;
  for(int i=0; i<porofluid_artificial_velocity_->MyLength();++i)
  {
    if(i%4==3)
    {
      int gid = poroelast_subproblem_->FluidField()->DofRowMap()->GID(i);
      values = porofluid_artificial_velocity_ ->Values();
      porofluid_artificial_source_->ReplaceGlobalValue(gid,0,values[i]);
    }
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellMigration::SetFieldDt()
{
  // set parameters for pre-simulation
  if(initialize_cell_ and Step()==0 and (IterationCounter())[0]==0)
  {
    if(myrank_==0)
    {
      std::cout<<"\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
      std::cout<<"++++++                 CELL INITIALIZATION STEP                   ++++++"<<std::endl;
      std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
      std::cout<<"dt="<<initial_timestep_<<"\n"<<std::endl;
    }
    SetDt(initial_timestep_);
    poroelast_subproblem_->SetDt(initial_timestep_);
    poroelast_subproblem_->PoroField()->SetDt(initial_timestep_);
    poroelast_subproblem_->FluidField()->SetDt(initial_timestep_);
    porostructure_->SetDt(initial_timestep_);
    poroelast_subproblem_->ScaTraField()->SetDt(initial_timestep_);
    cellstructure_->SetDt(initial_timestep_);
  }
  else if(initialize_cell_ and Step()==initialization_steps_)
  {
    if(myrank_==0)
    {
    std::cout<<"\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
    std::cout<<"++++++              CELL INITIALIZATION SUCCESSFUL                ++++++"<<std::endl;
    std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
    std::cout<<" RESET TIMESTEP TO dt="<<timestep_<<"\n"<<std::endl;
    }
    SetDt(timestep_);
    poroelast_subproblem_->SetDt(timestep_);
    poroelast_subproblem_->PoroField()->SetDt(timestep_);
    poroelast_subproblem_->FluidField()->SetDt(timestep_);
    porostructure_->SetDt(timestep_);
    poroelast_subproblem_->ScaTraField()->SetDt(timestep_);
    cellstructure_->SetDt(timestep_);
    initialize_cell_ = 0;

    double simpleECMInteractionConstant = ((((exchange_manager_->GetPointerToPenaltyTractionAtGPMap()->find(3))->second).at(0)).Norm2())/0.2;
    exchange_manager_->SetSimpleECMInteractionConstant( simpleECMInteractionConstant );
    std::cout<<"Set Slope of ECM Interaction Force to "<<simpleECMInteractionConstant<<std::endl;
  }

  return;
}


void IMMERSED::ImmersedPartitionedCellMigration::ReadRestart(int step)
{
  cellstructure_->ReadRestart(step);
  poroelast_subproblem_->ReadRestart(step);
  SetTimeStep(poroelast_subproblem_->PoroField()->Time(),step);

  return;
}
///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//void IMMERSED::ImmersedPartitionedCellMigration::DeleteForceAtLeadingEdge(Teuchos::RCP<Epetra_Vector> force)
//{
//
//  std::multimap<std::string,Teuchos::RCP<DRT::Condition> >::iterator condition_iterator;
//  //DRT::Condition* condition = immerseddis_->GetCondition("SurfaceNeumann");
//  std::multimap<std::string,Teuchos::RCP<DRT::Condition> > conditions = immerseddis_->GetAllConditions();
//
//  for (condition_iterator=conditions.begin(); condition_iterator!=conditions.end(); ++condition_iterator)
//  {
//    if(condition_iterator->first == (std::string)"SurfaceNeumann")
//    {
//      const std::vector<int>* conditioned_nodes = condition_iterator->second->Nodes();
//      int num_conditioned_nodes = conditioned_nodes->size();
//      for(int i=0; i<num_conditioned_nodes; ++i)
//      {
//        DRT::Node*  node = immerseddis_->gNode((*conditioned_nodes)[i]);
//        std::vector<int> dofs = immerseddis_->Dof(node);
//        for(int numdof=0;numdof<(int)dofs.size();numdof++)
//        {
//          int err = force-> ReplaceGlobalValue(dofs[numdof],0,0.0);
//          if(err != 0)
//            dserror("ReplaceGlobalValue returned %d for node ID",err,(*conditioned_nodes)[i]);
//        }
//      }
//    } // if surfaceneumann condition
//  } // loop over all conditions
//
//  deleted_force_leading_edge_ = true;
//  return;
//}
//
//
///*----------------------------------------------------------------------*/
///*----------------------------------------------------------------------*/
//void IMMERSED::ImmersedPartitionedCellMigration::PreventInterfacePullingTraction(Teuchos::RCP<Epetra_Vector> increment_vector,
//                                                                                 Teuchos::RCP<Epetra_Vector> total_vector)
//{
//
//  int mylength=increment_vector->MyLength();
//  const Epetra_BlockMap& map = increment_vector->Map();
//  int gid = -1234;
//  int numdofpernode = 3;
//  int numnodes = mylength/numdofpernode;
//
//  LINALG::Matrix<3,1> total_traction(true);
//  LINALG::Matrix<3,1> increment_traction(true);
//
//  double* increment_values = increment_vector->Values();
//  double* total_values = total_vector->Values();
//
//  double same_direction          = 1234.0;
//  double norm_increment_traction =-1234.0;
//  double norm_total_traction     =-1234.0;
//
//  // loop over all nodes
//  for(int node=0;node<numnodes;node++)
//  {
//
//    // loop over dofs
//    for(int dof=0;dof<numdofpernode;++dof)
//    {
//      increment_traction(dof)=increment_values[node*numdofpernode+dof];
//      total_traction(dof)=total_values[node*numdofpernode+dof];
//    } // dof loop
//
//
//    norm_increment_traction = increment_traction.Norm2();
//    norm_total_traction = total_traction.Norm2();
//
//    if(norm_total_traction>1.0e-09 and norm_increment_traction>1.0e-09)
//    {
//      same_direction = total_traction.Dot(increment_traction);
//
//      // only if increment is pointing in opposite direction of the total traction
//      if(same_direction < 0.0)
//      {
//        // only if the considered increment vector is additionally larger than the total traction
//        if(norm_increment_traction>norm_total_traction)
//        { std::cout<<"WARNING! PreventInterfacePullingTraction() had to zero interaction force and increment ... "<<std::endl;
//          // 2nd loop over dofs
//          for(int dof=0;dof<numdofpernode;++dof)
//          {
//            gid = map.GID(node*numdofpernode+dof);
//            total_vector->ReplaceGlobalValue(gid,0,0.0);
//            increment_vector->ReplaceGlobalValue(gid,0,0.0);
//          }
//        }
//      }
//
//    }
//  } // node loop
//
//  return;
//}
