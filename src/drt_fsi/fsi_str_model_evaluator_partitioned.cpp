/*-----------------------------------------------------------*/
/*! \file

\brief Model evaluator for structure part of partitioned fsi

\maintainer Nora Hagmeyer

\level 3

*/
/*-----------------------------------------------------------*/


#include "fsi_str_model_evaluator_partitioned.H"

#include "../drt_structure_new/str_model_evaluator.H"
#include "../drt_structure_new/str_dbc.H"
#include "../drt_structure_new/str_impl_generic.H"
#include "../drt_structure_new/str_timint_implicit.H"
#include "../drt_structure_new/str_nln_solver_generic.H"
#include "../drt_structure_new/str_timint_basedataglobalstate.H"

#include "../solver_nonlin_nox/nox_nln_group.H"
#include "../solver_nonlin_nox/nox_nln_aux.H"

#include "../linalg/linalg_utils.H"

#include "Epetra_Comm.h"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
STR::MODELEVALUATOR::PartitionedFSI::PartitionedFSI()
    : interface_force_np_ptr_(Teuchos::null), is_relaxationsolve(false)
{
  // empty
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedFSI::Setup()
{
  // fsi interface force at t_{n+1}
  interface_force_np_ptr_ = Teuchos::rcp(new Epetra_Vector(*GState().DofRowMap(), true));

  // set flag
  issetup_ = true;

  return;
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedFSI::SetupMultiMapExtractor()
{
  Int().ModelEval().SetupMultiMapExtractor();
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> STR::MODELEVALUATOR::PartitionedFSI::GetBlockDofRowMapPtr() const
{
  CheckInitSetup();
  return GState().DofRowMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::PartitionedFSI::GetCurrentSolutionPtr() const
{
  CheckInit();
  return GState().GetDisNp();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::PartitionedFSI::GetLastTimeStepSolutionPtr()
    const
{
  CheckInit();
  return GState().GetDisN();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::PartitionedFSI::AssembleForce(
    Epetra_Vector& f, const double& timefac_np) const
{
  LINALG::AssembleMyVector(1.0, f, -timefac_np, *interface_force_np_ptr_);
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::PartitionedFSI::UpdateStepState(const double& timefac_n)
{
  if (not is_relaxationsolve)  // standard case
  {
    // add the old time factor scaled contributions to the residual
    Teuchos::RCP<Epetra_Vector>& fstructold_ptr = GState().GetMutableFstructureOld();
    fstructold_ptr->Update(-timefac_n, *interface_force_np_ptr_, 1.0);
  }
  else
  {
    // do nothing
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> STR::MODELEVALUATOR::PartitionedFSI::SolveRelaxationLinear(
    Teuchos::RCP<ADAPTER::Structure> structure)
{
  // print to screen
  if (GState().DofRowMap()->Comm().MyPID() == 0)
    std::cout << "\n DO SRUCTURAL RELAXATION SOLVE ..." << std::endl;

  // cast adapter structure to implicit time integrator
  Teuchos::RCP<STR::TIMINT::Implicit> ti_impl =
      Teuchos::rcp_dynamic_cast<STR::TIMINT::Implicit>(structure, true);

  // get the nonlinear solver pointer
  STR::NLN::SOLVER::Generic& nlnsolver =
      const_cast<STR::NLN::SOLVER::Generic&>(*(ti_impl->GetNlnSolverPtr()));

  // get the solution group
  NOX::Abstract::Group& grp = nlnsolver.SolutionGroup();
  NOX::NLN::Group* grp_ptr = dynamic_cast<NOX::NLN::Group*>(&grp);
  if (grp_ptr == NULL) dserror("Dynamic cast failed!");

  // get nox parameter
  Teuchos::ParameterList& noxparams = ti_impl->DataSDyn().GetMutableNoxParams();

  // create new state vector
  Teuchos::RCP<NOX::Epetra::Vector> x_ptr = GState().CreateGlobalVector(
      DRT::UTILS::vec_init_last_time_step, ti_impl->ImplIntPtr()->ModelEvalPtr());
  // Set the solution vector in the nox group. This will reset all isValid
  // flags.
  grp.setX(*x_ptr);

  // ---------------------------------------------------------------------------
  // Compute F and jacobian
  // ---------------------------------------------------------------------------
  grp_ptr->computeJacobian();

  // overwrite F with boundary force
  interface_force_np_ptr_->Scale(-(ti_impl->TimIntParam()));
  ti_impl->DBCPtr()->ApplyDirichletToRhs(interface_force_np_ptr_);
  Teuchos::RCP<NOX::Epetra::Vector> nox_force =
      Teuchos::rcp(new NOX::Epetra::Vector(interface_force_np_ptr_));
  grp_ptr->setF(nox_force);

  // ---------------------------------------------------------------------------
  // Check if we are using a Newton direction
  // ---------------------------------------------------------------------------
  const std::string dir_str(NOX::NLN::AUX::GetDirectionMethodListName(noxparams));
  if (dir_str != "Newton")
    dserror(
        "The RelaxationSolve is currently only working for the direction-"
        "method \"Newton\".");

  // ---------------------------------------------------------------------------
  // (re)set the linear solver parameters
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList& p =
      noxparams.sublist("Direction").sublist("Newton").sublist("Linear Solver");
  p.set<int>("Number of Nonlinear Iterations", 0);
  p.set<int>("Current Time Step", GState().GetStepNp());
  p.set<double>("Wanted Tolerance", 1.0e-6);  //!< dummy

  // ---------------------------------------------------------------------------
  // solve the linear system of equations and update the current state
  // ---------------------------------------------------------------------------
  // compute the Newton direction
  grp_ptr->computeNewton(p);

  // get the increment from the previous solution step
  const NOX::Epetra::Vector& increment =
      dynamic_cast<const NOX::Epetra::Vector&>(grp_ptr->getNewton());

  // return the increment
  return Teuchos::rcpFromRef(increment.getEpetraVector());
}
