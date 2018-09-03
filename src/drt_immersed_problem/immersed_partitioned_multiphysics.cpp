/*!----------------------------------------------------------------------
\file immersed_partitioned_multiphysics.cpp

\brief base class for all multifield partitioned immersed algorithms

\level 2

<pre>
\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240
</pre>
*----------------------------------------------------------------------*/
#include "immersed_partitioned_multiphysics.H"

#include "immersed_partitioned_adhesion_traction.H"
#include "immersed_partitioned_flow_cell_interaction.H"

#include "../drt_ssi/ssi_partitioned_2wc.H"
#include "../drt_poroelast/poro_scatra_base.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_adapter/ad_str_multiphysicswrapper_cellmigration.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
IMMERSED::ImmersedPartitionedMultiphysics::ImmersedPartitionedMultiphysics(const Epetra_Comm& comm)
    : ImmersedPartitioned(comm),
      myrank_(comm.MyPID()),
      globalproblem_(NULL),
      immerseddis_(Teuchos::null),
      backgroundstructuredis_(Teuchos::null),
      cfi_module_(Teuchos::null),
      adh_module_(Teuchos::null),
      cellstructure_(Teuchos::null),
      cellscatra_subproblem_(Teuchos::null),
      poroscatra_subproblem_(Teuchos::null)
{
  // empty constructor
  return;
}  // ImmersedPartitionedMultiphysics


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int IMMERSED::ImmersedPartitionedMultiphysics::Init(const Teuchos::ParameterList& params)
{
  // reset the setup flag
  SetIsSetup(false);

  // get pointer to global problem
  globalproblem_ = DRT::Problem::Instance();

  // set pointer to discretization
  immerseddis_ = globalproblem_->GetDis("cell");
  backgroundfluiddis_ = globalproblem_->GetDis("porofluid");
  backgroundstructuredis_ = globalproblem_->GetDis("structure");


  // extract ptrs to modules from parameter list
  cfi_module_ = params.get<Teuchos::RCP<IMMERSED::ImmersedPartitionedFlowCellInteraction>>(
      "ImmersedPartitionedFlowCellInteraction");
  adh_module_ = params.get<Teuchos::RCP<IMMERSED::ImmersedPartitionedAdhesionTraction>>(
      "ImmersedPartitionedAdhesionTraction");

  if (cfi_module_ != Teuchos::null) cfi_module_->Init(params);
  if (adh_module_ != Teuchos::null) adh_module_->Init(params);

  // get pointer to cell structure
  Teuchos::RCP<ADAPTER::MultiphysicsStructureWrapperCellMigration> multiphysicswrapper =
      params.get<Teuchos::RCP<ADAPTER::MultiphysicsStructureWrapperCellMigration>>(
          "RCPToCellStructure");
  if (multiphysicswrapper == Teuchos::null)
    dserror("no pointer to MultiphysicsStructureWrapperCellMigration provided");
  cellstructure_ = multiphysicswrapper->GetFSIStructureWrapperPtr();

  // get pointer structure-scatra interaction (ssi) subproblem
  cellscatra_subproblem_ = params.get<Teuchos::RCP<SSI::SSI_Part2WC>>("RCPToCellScatra");

  // get pointer poroelast-scatra interaction (ssi) subproblem
  poroscatra_subproblem_ = params.get<Teuchos::RCP<POROELAST::PoroScatraBase>>("RCPToPoroScatra");

  // set isinit_ flag true
  SetIsInit(true);

  return 0;
}  // Init


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedMultiphysics::Setup()
{
  // make sure Init(...) was called first
  CheckIsInit();

  if (cfi_module_ != Teuchos::null) cfi_module_->Setup();
  if (adh_module_ != Teuchos::null) adh_module_->Setup();

  // set flag issetup true
  SetIsSetup(true);
}  // Setup


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedMultiphysics::Solve()
{
  while (NotFinished())
  {
    PrepareTimeStep();

    TimeStep();

    PrepareOutput();

    Update();

    Output();
  }

  return;
}  // Solve


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedMultiphysics::TimeStep()
{
  CFIOperator();
  CEIOperator();

  return;
}  // TimeStep


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedMultiphysics::CFIOperator()
{
  // declare epetra map
  Teuchos::RCP<Epetra_Map> dofmap_cell = Teuchos::null;
  Teuchos::RCP<Epetra_Map> dofmap_ecm = Teuchos::null;

  // apply Dirichlet to adhesion nodes
  ApplyDirichlet(
      cellstructure_, immerseddis_, "CellFocalAdhesion", dofmap_cell, 3, cellstructure_->Dispnp());

  // set state in nox and gobal state object
  cellstructure_->SetState(cellstructure_->WriteAccessDispnp());

  // apply Dirichlet to ECM structure
  ApplyDirichlet(poroscatra_subproblem_->PoroField()->StructureField(), backgroundstructuredis_,
      "PoroCoupling", dofmap_ecm, 3,
      poroscatra_subproblem_->PoroField()->StructureField()->Dispnp());

  // rebuild the monolithic poro dbc map
  poroscatra_subproblem_->PoroField()->BuildCombinedDBCMap();
  // solve cfi module for reaction forces at adhesion nodes
  Teuchos::RCP<Epetra_Vector> cell_bdry_traction = cfi_module_->DoStep(cfi_module_, Teuchos::null);

  // remove Dirichlet condition from adhesion nodes
  RemoveDirichlet(dofmap_cell, cellstructure_);
  RemoveDirichlet(dofmap_ecm, poroscatra_subproblem_->PoroField()->StructureField());

  return;
}  // CFIOperator


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedMultiphysics::CEIOperator()
{
  // declare epetra map
  Teuchos::RCP<Epetra_Map> dofmap_ecm = Teuchos::null;

  //   // apply Dirichlet to interstitial flow
  //   ApplyDirichletToFluid(
  //       poroscatra_subproblem_->FluidField(),
  //       backgroundfluiddis_,
  //       "PoroCoupling",
  //       dofmap_ecm,
  //       4,
  //       poroscatra_subproblem_->FluidField()->Velnp());

  //   const double dt = globalproblem_->CellMigrationParams().get<double>("TIMESTEP");
  //
  //   disp_->Scale(0.0);
  //   disp_->Update(1.0,*poroscatra_subproblem_->PoroField()->StructureField()->Dispn(),1.0);
  //   disp_->Update(dt,
  //       *poroscatra_subproblem_->PoroField()->FluidToStructureField(poroscatra_subproblem_->FluidField()->Velnp()),
  //       1.0);
  //
  //   ApplyDirichletToArtificialECM(
  //       poroscatra_subproblem_->PoroField()->StructureField(),
  //       backgroundstructuredis_,
  //       "PoroCoupling",
  //       dofmap_ecm,
  //       3,
  //       disp_);

  // rebuild the monolithic poro dbc map
  poroscatra_subproblem_->PoroField()->BuildCombinedDBCMap();
  // solve adhesion module
  adh_module_->DoStep(adh_module_, Teuchos::null);

  //   RemoveDirichletFromFluid(dofmap_ecm,poroscatra_subproblem_->FluidField());
  //   RemoveDirichlet(dofmap_ecm,poroscatra_subproblem_->PoroField()->StructureField());

  return;
}  // CEIOperator


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedMultiphysics::PrepareTimeStep()
{
  IncrementTimeAndStep();

  if (myrank_ == 0)
  {
    std::cout << "+------------------------------------+" << std::endl;
    std::cout << "+ Prepare Multiphysics Time Step ... |" << std::endl;
    std::cout << "+------------------------------------+" << std::endl;
    std::cout << "Step=" << Step() << " Time=" << Time() << " dt=" << Dt() << std::endl;
  }
  cellscatra_subproblem_->PrepareTimeStep(false);
  poroscatra_subproblem_->PrepareTimeStep(false);

  return;
}  // PrepareTimeStep

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedMultiphysics::PrepareOutput()
{
  cellscatra_subproblem_->StructureField()->PrepareOutput();
  poroscatra_subproblem_->PrepareOutput();

  return;
}  // PrepareOutput

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedMultiphysics::Update()
{
  cellscatra_subproblem_->StructureField()->Update();
  cellscatra_subproblem_->ScaTraField()->ScaTraField()->Update();
  poroscatra_subproblem_->Update();

  return;
}  // Update


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedMultiphysics::Output()
{
  cellscatra_subproblem_->StructureField()->Output();
  cellscatra_subproblem_->ScaTraField()->ScaTraField()->Output();
  poroscatra_subproblem_->Output();

  cfi_module_->WriteDrag();

  return;
}  // Output



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedMultiphysics::ReadRestart(int step) { return; }  // ReadRestart
