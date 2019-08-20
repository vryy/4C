/*----------------------------------------------------------------------*/
/*! \file

\brief specialization of ssi2wc, including adhesion dynamics

\level 2

\maintainer Jonas Eichinger
 *----------------------------------------------------------------------*/


#include "ssi_partitioned_2wc_adhesiondynamics.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_base.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/ad_str_ssiwrapper.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io_control.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"



/*----------------------------------------------------------------------*
 | constructor                                                          |
 *----------------------------------------------------------------------*/
SSI::SSI_Part2WC_ADHESIONDYNAMICS::SSI_Part2WC_ADHESIONDYNAMICS(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSI_Part2WC(comm, globaltimeparams),
      cell_adhesion_forces_(Teuchos::null),
      cell_adhesion_fixpoints_(Teuchos::null),
      exchange_manager_(DRT::ImmersedFieldExchangeManager::Instance())
{
  // empty
}


/*----------------------------------------------------------------------*
 | Initialize this object                                   rauch 08/16 |
 *----------------------------------------------------------------------*/
int SSI::SSI_Part2WC_ADHESIONDYNAMICS::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams, const std::string struct_disname,
    const std::string scatra_disname, bool isAle)
{
  int returnvar = 0;

  // call setup of base class
  returnvar = SSI::SSI_Part2WC::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);

  // check if scatra in cell is set up with ale description
  if (not ScaTraField()->ScaTraField()->IsALE())
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
  BuildConditionDofRowMap(StructureField()->Discretization()->GetCondition("CellFocalAdhesion"),
      StructureField()->Discretization(), conditiondofrowmap_);

  // build condition dof map with row layout
  BuildConditionDofColMap(StructureField()->Discretization()->GetCondition("CellFocalAdhesion"),
      StructureField()->Discretization(), conditiondofcolmap_);

  // create traction vector
  surface_traction_ = LINALG::CreateVector(*conditiondofrowmap_, true);
  surface_traction_col_ = LINALG::CreateVector(*conditiondofcolmap_, true);

  // give traction vector pointer to ExchangeManager
  exchange_manager_->SetPointerSurfaceTraction(surface_traction_col_);
}


/*------------------------------------------------------------------------*
 | BuildConditionDofRowMap                                    rauch 03/16 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ADHESIONDYNAMICS::BuildConditionDofRowMap(const DRT::Condition* condition,
    const Teuchos::RCP<const DRT::Discretization>& dis, Teuchos::RCP<Epetra_Map>& condnodemap)
{
  const std::vector<int>* conditionednodes = condition->Nodes();

  std::vector<int> mydirichdofs(0);

  for (int i = 0; i < (int)conditionednodes->size(); ++i)
  {
    DRT::Node* currnode = dis->gNode(conditionednodes->at(i));
    std::vector<int> dofs = dis->Dof(0, currnode);

    for (int dim = 0; dim < 3; ++dim)
    {
      if (dis->DofRowMap()->LID(dofs[dim]) != -1) mydirichdofs.push_back(dofs[dim]);
    }
  }

  int nummydirichvals = mydirichdofs.size();
  condnodemap =
      Teuchos::rcp(new Epetra_Map(-1, nummydirichvals, &(mydirichdofs[0]), 0, dis->Comm()));

  return;
}


/*------------------------------------------------------------------------*
 | BuildConditionDofColMap                                    rauch 03/16 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ADHESIONDYNAMICS::BuildConditionDofColMap(const DRT::Condition* condition,
    const Teuchos::RCP<const DRT::Discretization>& dis, Teuchos::RCP<Epetra_Map>& condnodemap)
{
  const std::vector<int>* conditionednodes = condition->Nodes();

  std::vector<int> mydirichdofs(0);

  for (int i = 0; i < (int)conditionednodes->size(); ++i)
  {
    DRT::Node* currnode = dis->gNode(conditionednodes->at(i));
    std::vector<int> dofs = dis->Dof(0, currnode);

    for (int dim = 0; dim < 3; ++dim)
    {
      if (dis->DofColMap()->LID(dofs[dim]) != -1) mydirichdofs.push_back(dofs[dim]);
    }
  }

  int nummydirichvals = mydirichdofs.size();
  condnodemap =
      Teuchos::rcp(new Epetra_Map(-1, nummydirichvals, &(mydirichdofs[0]), 0, dis->Comm()));

  return;
}


/*------------------------------------------------------------------------*
 | update the current states in every iteration               rauch 05/16 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ADHESIONDYNAMICS::IterUpdateStates()
{
  // perform the update from the base class
  SSI::SSI_Part2WC::IterUpdateStates();

  // clear surface traction vector
  surface_traction_->Scale(0.0);

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
  surface_traction_->Scale(0.0);

  return;
}


/*------------------------------------------------------------------------*
 | Pre operator for second field operator                     rauch 06/17 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ADHESIONDYNAMICS::PreOperator2()
{
  if (iter_ != 1 and use_old_structure_)
  {
    // NOTE: the predictor is NOT called in here. Just the screen output is not correct.
    // we only get norm of the evaluation of the structure problem
    structure_->PreparePartitionStep();
  }

  return;
}


/*------------------------------------------------------------------------*
 | Post operator for first field operator                     rauch 06/17 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ADHESIONDYNAMICS::PostOperator1()
{
  Teuchos::RCP<Epetra_Vector> phinp = scatra_->ScaTraField()->Phinp();
  double norm = 0.0;
  phinp->Norm2(&norm);


  if (Comm().MyPID() == 0)
  {
    std::cout << "L2-Norm of Concentration Vector: " << std::setprecision(11) << norm << std::endl;
    std::cout << "---------------" << std::endl;
  }

  return;
}


/*------------------------------------------------------------------------*
 | Post operator for second field operator                    rauch 06/17 |
 *------------------------------------------------------------------------*/
void SSI::SSI_Part2WC_ADHESIONDYNAMICS::PostOperator2()
{
  // print norm of cell adhesion force vector
  // first we have to get the adhesion vector. this can't be done in Setup(), as the pointer is not
  // yet set there...
  if (cell_adhesion_forces_ == Teuchos::null)
    cell_adhesion_forces_ = exchange_manager_->GetPointerCellAdhesionForce();
  double fadhnorm = -1234.0;
  cell_adhesion_forces_->Norm2(&fadhnorm);
  if (Comm().MyPID() == 0)
  {
    std::cout << std::endl
              << "L2-Norm of Cell Adhesion Force Vector: " << std::setprecision(11) << fadhnorm
              << std::endl;
    std::cout << "---------------" << std::endl;
  }

  return;
}
