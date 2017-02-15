  /*!----------------------------------------------------------------------
\file immersed_partitioned_cellcontraction.cpp

\brief partitioned immersed contraction algorithm

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

\level 3

*----------------------------------------------------------------------*/
#include "immersed_partitioned_cellcontraction.H"

#include "ssi_partitioned_2wc_biochemomechano.H"

#include "../drt_poroelast/poro_scatra_base.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_adapter/ad_ale_fluid.H"
#include "../drt_adapter/ad_str_fsiwrapper_immersed.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_str_multiphysicswrapper_cellmigration.H"

IMMERSED::ImmersedPartitionedCellContraction::ImmersedPartitionedCellContraction(
    const Teuchos::ParameterList& params,
    const Epetra_Comm& comm)
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
  Teuchos::RCP<ADAPTER::MultiphysicsStructureWrapperCellMigration> multiphysicswrapper =
      params.get<Teuchos::RCP<ADAPTER::MultiphysicsStructureWrapperCellMigration> >("RCPToCellStructure");

  if(multiphysicswrapper == Teuchos::null)
    dserror("no pointer to MultiphysicsStructureWrapperCellMigration provided");

  cellstructure_ = multiphysicswrapper->GetFSIStructureWrapperPtr();

  // get pointer poroelast-scatra interaction subproblem
  poroscatra_subproblem_ = params.get<Teuchos::RCP<POROELAST::PoroScatraBase> >("RCPToPoroScatra");

  // get pointer structure-scatra interaction (ssi) subproblem
  Teuchos::RCP<SSI::SSI_Part2WC> cellscatra_subproblem = params.get<Teuchos::RCP<SSI::SSI_Part2WC> >("RCPToCellScatra");
  cellscatra_subproblem_ = Teuchos::rcp_dynamic_cast<SSI::SSI_Part2WC_BIOCHEMOMECHANO>(cellscatra_subproblem);

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
  if(cellscatra_subproblem_==Teuchos::null)
    dserror("no pointer to cellscatra_subproblem_ provided !");

  // initialize important variables for parallel simulations
  myrank_  = comm.MyPID();
  numproc_ = comm.NumProc();

  // get pointer to global problem
  globalproblem_ = DRT::Problem::Instance();

  // get pointer to discretizations
  backgroundfluiddis_     = globalproblem_->GetDis("porofluid");
  backgroundstructuredis_ = globalproblem_->GetDis("structure");
  immerseddis_            = globalproblem_->GetDis("cell");
  scatradis_              = globalproblem_->GetDis("scatra");

  // get coupling variable
  displacementcoupling_ = globalproblem_->ImmersedMethodParams().sublist("PARTITIONED SOLVER").get<std::string>("COUPVARIABLE") == "Displacement";
  if(displacementcoupling_ and myrank_==0)
    std::cout<<"\n Coupling variable for partitioned protrusion formation:  Displacements "<<std::endl;
  else if (!displacementcoupling_ and myrank_==0)
    std::cout<<"\n Coupling variable for partitioned protrusion formation :  Force "<<std::endl;

  // construct immersed exchange manager. singleton class that makes immersed variables comfortably accessible from everywhere in the code
  exchange_manager_ = DRT::ImmersedFieldExchangeManager::Instance();

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int IMMERSED::ImmersedPartitionedCellContraction::Init(const Teuchos::ParameterList& params)
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
void IMMERSED::ImmersedPartitionedCellContraction::Setup()
{
  // make sure Init(...) was called first
  CheckIsInit();

  // get parameters for nox
  const Teuchos::ParameterList& immerseddyn =
      globalproblem_->ImmersedMethodParams().sublist("PARTITIONED SOLVER");
  SetDefaultParameters(immerseddyn,NOXParameterList());
  //noxparameterlist_.print();

  // set flag issetup true
  SetIsSetup(true);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellContraction::ReadRestart(int step)
{
  cellstructure_->ReadRestart(step);
  poroscatra_subproblem_->ReadRestart(step);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedCellContraction::InitialGuess()
{
  if(displacementcoupling_)
    return cellstructure_->PredictImmersedInterfaceDispnp();
  else // FORCE COUPLING
  {
    dserror("Force Coupling for Cell Contraction is not implemented, yet.");
    return Teuchos::null;
  }

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellContraction::PrepareTimeStep()
{

  IncrementTimeAndStep();

  PrintHeader();

  cellscatra_subproblem_->PrepareTimeStep(false);

  poroscatra_subproblem_->PrepareTimeStep(false);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellContraction::CouplingOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  // DISPLACEMENT COUPLING
  if (displacementcoupling_)
  {
    const Teuchos::RCP<Epetra_Vector> idispn = Teuchos::rcp(new Epetra_Vector(x));

    ////////////////////
    // CALL ImmersedOp
    ////////////////////
    PrepareImmersedOp();
    Teuchos::RCP<Epetra_Vector> idispnp =
    ImmersedOp(Teuchos::null, fillFlag);

    ////////////////////
    // CALL BackgroundOp
    ////////////////////
    PrepareBackgroundOp();
    BackgroundOp(Teuchos::null, fillFlag);

    int err = F.Update(1.0, *idispnp, -1.0, *idispn, 0.0);
    if(err != 0)
      dserror("Vector update of Coupling-residual returned err=%d",err);

  }
  // FORCE COUPLING
  else if(!displacementcoupling_)
  {
   dserror("Force Coupling for Protrusion Formation is not implemented, yet.");

  } // displacement / force coupling

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellContraction::BackgroundOp(Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values,
                                                              const FillType fillFlag)
{
  IMMERSED::ImmersedPartitioned::BackgroundOp(backgrd_dirichlet_values,fillFlag);

  if (fillFlag==User)
  {
    dserror("fillFlag == User : not yet implemented");
  }
  else
  {
    if(myrank_ == 0)
      std::cout<<"BackgroundOp is empty. So far, this is only a one way coupled Problem.\n"<<std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
IMMERSED::ImmersedPartitionedCellContraction::ImmersedOp(
    Teuchos::RCP<Epetra_Vector> bdry_traction,
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

     // solve cell
    if (Comm().MyPID()==0)
    {
      std::cout<<"\n****************************************\n          PROTRUSION FORMATION\n****************************************\n";
    }
    cellscatra_subproblem_->OuterLoop();

    return cellstructure_->ExtractImmersedInterfaceDispnp();
  }

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellContraction::PrepareBackgroundOp()
{
  // do nothing so far

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellContraction::PrepareImmersedOp()
{
 // do nothing so far

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellContraction::PrepareOutput()
{
  poroscatra_subproblem_->PrepareOutput();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellContraction::Output()
{
  poroscatra_subproblem_->Output();
  cellscatra_subproblem_->UpdateAndOutput();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitionedCellContraction::Update()
{
  poroscatra_subproblem_->Update();

  return;
}
