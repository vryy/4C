/*-----------------------------------------------------------*/
/*!
\file nox_nln_problem_factory.cpp

\maintainer Michael Hiermeier

\date Jul 23, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_problem_factory.H"    // class definition
#include "nox_nln_constraint_interface_required.H"
#include "nox_nln_globaldata.H"

#include <Epetra_Vector.h>

#include "../linalg/linalg_sparseoperator.H"

// all the different NoxProblems
#include "nox_nln_problem.H"
#include "nox_nln_constraint_problem.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::Problem::Factory::Factory()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::NoxProblem> NOX::NLN::Problem::Factory::BuildNoxProblem(
    Teuchos::RCP<Epetra_Vector>& x,
    Teuchos::RCP<LINALG::SparseOperator> A,
    Teuchos::RCP<NOX::NLN::GlobalData>& nlnGlobalData)
{
  Teuchos::RCP<NOX::NLN::NoxProblem> noxProblem;

  // -------------------------------------
  // constrained optimization problems
  // -------------------------------------
  if (nlnGlobalData->GetIsConstrained())
    noxProblem = Teuchos::rcp(new NOX::NLN::CONSTRAINT::NoxProblem(nlnGlobalData,x,A));
  // -------------------------------------
  // unconstrained optimization problems
  // -------------------------------------
  else
    noxProblem = Teuchos::rcp(new NOX::NLN::NoxProblem(nlnGlobalData,x,A));

  return noxProblem;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::NLN::NoxProblem> NOX::NLN::Problem::BuildNoxProblem(
        Teuchos::RCP<Epetra_Vector>& x,
        Teuchos::RCP<LINALG::SparseOperator> A,
        Teuchos::RCP<NOX::NLN::GlobalData>& nlnGlobalData)
{
  Factory factory;
  return factory.BuildNoxProblem(x,A,nlnGlobalData);
}

