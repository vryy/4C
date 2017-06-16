/*!----------------------------------------------------------------------
\file ssi_partitioned_2wc_adhesiondynamics.cpp

\brief specialization of ssi2wc, including adhesion dynamics

\level 3

\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15240

 *----------------------------------------------------------------------*/


#include "ssi_partitioned_2wc_adhesiondynamics.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_base.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../linalg/linalg_utils.H"

#include "../drt_io/io_control.H"
#include "../linalg/linalg_solver.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 | constructor                                                          |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_ADHESIONDYNAMICS::SSI_Part2WC_ADHESIONDYNAMICS(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams)
: SSI_Part2WC(comm, globaltimeparams),
  exchange_manager_(DRT::ImmersedFieldExchangeManager::Instance())
{
  // empty
}


/*----------------------------------------------------------------------*
 | Setup this object                                        rauch 08/16 |
 *----------------------------------------------------------------------*/
int SSI::SSI_Part2WC_ADHESIONDYNAMICS::Init(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams,
    const std::string struct_disname,
    const std::string scatra_disname,
    bool isAle)
{
  int returnvar=0;

  // call setup of base class
  returnvar =
      SSI::SSI_Part2WC::Init(comm,globaltimeparams,scatraparams,structparams,struct_disname,scatra_disname,isAle);

  // check if scatra in cell is set up with ale description
  if(not ScaTraField()->ScaTraField()->IsALE())
    dserror("We need an ALE description for the cell-scatra field!");

  return returnvar;
}


/*------------------------------------------------------------------------*
 | Setup this object                                          rauch 12/16 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ADHESIONDYNAMICS::Setup()
{
  // call setup of base class
  SSI::SSI_Part2WC::Setup();

  // build condition dof map with row layout
  BuildConditionDofRowMap(
      StructureField()->Discretization()->GetCondition("CellFocalAdhesion"),
      StructureField()->Discretization(),
      conditiondofrowmap_);

  // build condition dof map with row layout
  BuildConditionDofColMap(
      StructureField()->Discretization()->GetCondition("CellFocalAdhesion"),
      StructureField()->Discretization(),
      conditiondofcolmap_);

  // create traction vector
  surface_traction_ = LINALG::CreateVector(*conditiondofrowmap_,true);
  surface_traction_col_ = LINALG::CreateVector(*conditiondofcolmap_,true);

  // give traction vector pointer to ExchangeManager
  exchange_manager_->SetPointerSurfaceTraction(surface_traction_col_);
}


/*------------------------------------------------------------------------*
 | BuildConditionDofRowMap                                    rauch 03/16 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ADHESIONDYNAMICS::BuildConditionDofRowMap(
    const DRT::Condition* condition,
    const Teuchos::RCP<const DRT::Discretization>& dis,
    Teuchos::RCP<Epetra_Map>& condnodemap)
{
  const std::vector<int>* conditionednodes = condition->Nodes();

  std::vector<int> mydirichdofs(0);

  for(int i=0; i<(int)conditionednodes->size(); ++i)
  {
    DRT::Node* currnode = dis->gNode(conditionednodes->at(i));
    std::vector<int> dofs = dis->Dof(0,currnode);

    for (int dim=0;dim<3;++dim)
    {
      if(dis->DofRowMap()->LID(dofs[dim]) != -1)
        mydirichdofs.push_back(dofs[dim]);
    }
  }

  int nummydirichvals = mydirichdofs.size();
  condnodemap = Teuchos::rcp( new Epetra_Map(-1,nummydirichvals,&(mydirichdofs[0]),0,dis->Comm()) );

  return;
}


/*------------------------------------------------------------------------*
 | BuildConditionDofColMap                                    rauch 03/16 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ADHESIONDYNAMICS::BuildConditionDofColMap(
    const DRT::Condition* condition,
    const Teuchos::RCP<const DRT::Discretization>& dis,
    Teuchos::RCP<Epetra_Map>& condnodemap)
{
  const std::vector<int>* conditionednodes = condition->Nodes();

  std::vector<int> mydirichdofs(0);

  for(int i=0; i<(int)conditionednodes->size(); ++i)
  {
    DRT::Node* currnode = dis->gNode(conditionednodes->at(i));
    std::vector<int> dofs = dis->Dof(0,currnode);

    for (int dim=0;dim<3;++dim)
    {
      if(dis->DofColMap()->LID(dofs[dim]) != -1)
        mydirichdofs.push_back(dofs[dim]);
    }
  }

  int nummydirichvals = mydirichdofs.size();
  condnodemap = Teuchos::rcp( new Epetra_Map(-1,nummydirichvals,&(mydirichdofs[0]),0,dis->Comm()) );

  return;
}


/*----------------------------------------------------------------------*
 | Solve Scatra field                                       rauch 12/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ADHESIONDYNAMICS::DoScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
    << "\n***********************\n  TRANSPORT SOLVER \n***********************\n";
  }

  // set displacement state
  StructureField()->Discretization()->SetState("displacement", StructureField()->Dispnp());
  // evaluate bond traction
  EvaluateBondTraction();

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->Solve();

  // set structure-based scalar transport values
  return SetScatraSolution(scatra_->ScaTraField()->Phinp());
}


/*----------------------------------------------------------------------*
 | Calculate traction at the adhesion surface               rauch 12/16 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ADHESIONDYNAMICS::EvaluateBondTraction()
{

  Teuchos::ParameterList params;

  // add action for bond traction evaluation
  params.set<std::string>("action","calc_cell_nodal_bond_traction");

  // least squares system-matrix
  Teuchos::RCP<LINALG::SparseOperator> leastsquares_matrix =
      Teuchos::rcp(new LINALG::SparseMatrix(*conditiondofrowmap_,18,false,true));
  // least squares right-hand side
  Teuchos::RCP<Epetra_Vector> leastsquares_rhs = LINALG::CreateVector(*conditiondofrowmap_,true);

  // evaluate least squares traction
  StructureField()->Discretization()->EvaluateCondition(params,
      leastsquares_matrix,
      Teuchos::null,
      leastsquares_rhs,
      Teuchos::null,
      Teuchos::null,
      "CellFocalAdhesion",
      -1);

  /////////////////////////////////////////////////////
  // Calc nodal traction from gauss point traction
  // with least squares method, i.e.:
  // Minimize sum|N * d - d_gp|²
  // And solve the resulting system of equations for
  // the leastsquares_traction x :
  // b = A * x --> x = A^1 * b
  // with
  // b = leastsquares_rhs from Evaluate Condition
  // A = leastsquares_matrix from Evaluate Condition
  ///////////////////////////////////////////////////
  // SOLVE LEAST SQUARES SYSTEM
  // solver setup
  Teuchos::ParameterList param_solve = DRT::Problem::Instance()->UMFPACKSolverParams();
  bool refactor = true;
  bool reset = true;
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(param_solve, Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
  StructureField()->Discretization()->ComputeNullSpaceIfNecessary(solver->Params());

  // complete matrix
  leastsquares_matrix->Complete();

  // solve for least squares optimal nodal values
  solver->Solve(leastsquares_matrix->EpetraOperator(), surface_traction_, leastsquares_rhs, refactor, reset);

  surface_traction_col_ -> Scale(0.0);

  // export to column map
  LINALG::Export(*surface_traction_,*surface_traction_col_);
}


/*------------------------------------------------------------------------*
 | update the current states in every iteration               rauch 05/16 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ADHESIONDYNAMICS::IterUpdateStates()
{
  // perform the update from the base class
  SSI::SSI_Part2WC::IterUpdateStates();

  // clear surface traction vector
  surface_traction_->PutScalar(0.0);

  return;
}


/*------------------------------------------------------------------------*
 | Update and Output                                          rauch 12/16 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ADHESIONDYNAMICS::UpdateAndOutput()
{
  // call update and output from base class
  SSI::SSI_Part2WC::UpdateAndOutput();

  // clear surface traction vector
  surface_traction_->PutScalar(0.0);

  return;
}
