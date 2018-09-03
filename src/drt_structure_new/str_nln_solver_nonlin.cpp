/*-----------------------------------------------------------*/
/*!
\file str_nln_solver_nonlin.cpp

\maintainer Matthias Mayr

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_nln_solver_nonlin.H"
#include "str_timint_implicit.H"

#include "nox_nln_str_linearsystem.H"

#include "../solver_nonlin/nln_problem.H"
#include "../solver_nonlin/nln_operator_base.H"
#include "../solver_nonlin/nln_operator_factory.H"

#include "../drt_inpar/inpar_structure.H"

#include "../linalg/linalg_sparseoperator.H"
#include "../linalg/linalg_solver.H"

#include <Epetra_Vector.h>

#include <Teuchos_ParameterList.hpp>

#include <NOX_Utils.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_Group.H>
#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Epetra_LinearSystem.H>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::NLN::SOLVER::Nonlin::Nonlin()
    : nlnproblem_(Teuchos::null),
      nlnoperator_(Teuchos::null),
      noxsoln_(Teuchos::null),
      noxutils_(Teuchos::null),
      params_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::Nonlin::Setup()
{
  CheckInit();

  // ---------------------------------------------------------------------------
  // Create / read parameter list for configuration of nonlinear solver
  // ---------------------------------------------------------------------------
  params_ = NLNSOL::UTILS::CreateParamListFromXML();

  // ---------------------------------------------------------------------------
  // Create the nox utils object
  // ---------------------------------------------------------------------------
  Teuchos::ParameterList printparams = *(NoxCreatePrintParameters(false));
  noxutils_ = Teuchos::rcp(new NOX::Utils(printparams));

  // set flag
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::NLN::SOLVER::Nonlin::Reset()
{
  CheckInitSetup();

  // do a hard reset at the beginning to be on the safe side
  nlnproblem_ = Teuchos::null;
  nlnoperator_ = Teuchos::null;
  GroupPtr() = Teuchos::null;
  noxsoln_ = Teuchos::null;

  // FixMe add the solution vector, right-hand-side vector and the complete
  // jacobian here!
  Teuchos::RCP<Epetra_Vector> soln = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null;
  Teuchos::RCP<LINALG::SparseOperator> jac = Teuchos::null;

  // ---------------------------------------------------------------------------
  // Create NOX group
  // ---------------------------------------------------------------------------
  // create initial guess vector of predictor result as CreateCopy to avoid
  // direct access to disn_
  noxsoln_ = Teuchos::rcp(new NOX::Epetra::Vector(soln, NOX::Epetra::Vector::CreateCopy));

  // create NOX linear system to provide access to Jacobian
  Teuchos::RCP<NOX::Epetra::LinearSystem> linsys = NoxCreateLinearSystem(*params_, *noxsoln_);

  /* use NOX::STR::Group to enable access to time integration
   * Note: NOX::Epetra::Group would be sufficient. */
  const Teuchos::RCP<NOX::Epetra::Interface::Required> ireq = NoxInterfacePtr();
  GroupPtr() = Teuchos::rcp(
      new NOX::Epetra::Group(params_->sublist("Nonlinear Problem"), ireq, *noxsoln_, linsys));

  // ---------------------------------------------------------------------------
  // Create interface to nonlinear problem
  // ---------------------------------------------------------------------------
  nlnproblem_ = Teuchos::rcp(new NLNSOL::NlnProblem());
  nlnproblem_->Init(
      DataGlobalState().GetComm(), params_->sublist("Nonlinear Problem"), *GroupPtr(), jac);
  nlnproblem_->Setup();

  /* Evaluate once more to guarantee valid quantities inside of the
   * NOX::STR::Group() */
  nlnproblem_->ComputeF(*soln, *rhs);
  nlnproblem_->ComputeJacobian();

  // ---------------------------------------------------------------------------
  // Create the nonlinear operator to solve the nonlinear problem
  // ---------------------------------------------------------------------------
  // use factory to create the nonlinear operator
  NLNSOL::NlnOperatorFactory opfactory;
  nlnoperator_ = opfactory.Create(params_->sublist("Nonlinear Operator"));

  // setup
  nlnoperator_->Init(
      DataGlobalState().GetComm(), params_->sublist("Nonlinear Operator"), nlnproblem_);
  nlnoperator_->Setup();

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum INPAR::STR::ConvergenceStatus STR::NLN::SOLVER::Nonlin::Solve()
{
  CheckInitSetup();
  // FixMe insert current right hand side
  Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null;
  int lnonlin_error = nlnoperator_->ApplyInverse(*rhs, noxsoln_->getEpetraVector());

  // Since it is possible that the nonlinear solution fails only on some procs
  // we need to communicate the error.
  int gnonlin_error = 0;
  DataGlobalState().GetComm().MaxAll(&lnonlin_error, &gnonlin_error, 1);

  INPAR::STR::ConvergenceStatus status = static_cast<INPAR::STR::ConvergenceStatus>(gnonlin_error);

  return status;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem> STR::NLN::SOLVER::Nonlin::NoxCreateLinearSystem(
    Teuchos::ParameterList& nlParams, NOX::Epetra::Vector& noxsoln)
{
  CheckInitSetup();

  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys = Teuchos::null;

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  Teuchos::RCP<NOX::Epetra::Interface::Required> ireq = NoxInterfacePtr();
  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> ijac = NoxInterfacePtr();

  // FixMe insert jacobian
  const Teuchos::RCP<LINALG::SparseOperator> jac = Teuchos::null;

  /* Support of only one linear solver for pure structural problems.
   * See the STR::NLN::SOLVER::Nox::ConvertModelType2SolType routine if you
   * consider to extend the functionality.                  hiermeier 10/2015
   */
  std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>> linsolver;
  linsolver[NOX::NLN::sol_structure] = DataSDyn().GetLinSolvers().at(INPAR::STR::model_structure);

  // call constructor without preconditioner and scaling object
  linSys = Teuchos::rcp(
      new NOX::NLN::STR::LinearSystem(printParams, lsParams, linsolver, ireq, ijac, jac, noxsoln));

  return linSys;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Teuchos::ParameterList> STR::NLN::SOLVER::Nonlin::NoxCreatePrintParameters(
    const bool verbose) const
{
  CheckInit();
  // Set the printing parameters in the "Printing" sublist
  Teuchos::RCP<Teuchos::ParameterList> printParams = Teuchos::rcp(new Teuchos::ParameterList());
  printParams->set("MyPID", DataGlobalState().GetMyRank());
  printParams->set("Output Precision", 6);
  printParams->set("Output Processor", 0);
  if (verbose)
  {
    printParams->set("Output Information",
        NOX::Utils::OuterIteration + NOX::Utils::OuterIterationStatusTest +
            NOX::Utils::InnerIteration + NOX::Utils::LinearSolverDetails + NOX::Utils::Parameters +
            NOX::Utils::Details + NOX::Utils::Warning + NOX::Utils::Debug +
            NOX::Utils::TestDetails + NOX::Utils::Error);
  }
  else if (ImplicitTimInt().GetDataIO().GetPrintIntermediateIterations())
  {
    printParams->set("Output Information",
        NOX::Utils::Error + NOX::Utils::OuterIterationStatusTest + NOX::Utils::TestDetails);
  }
  else
  {
    printParams->set("Output Information", NOX::Utils::Error + NOX::Utils::TestDetails);
  }

  // deliver liver
  return printParams;
}
