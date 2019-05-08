/*----------------------------------------------------------------------------*/
/**
\file xcontact_levelset_reinit_elliptic.cpp

\brief xcontact level-set elliptical reinitialization algorithm

\maintainer Matthias Mayr

\date Dec 1, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/


#include "xcontact_levelset_reinit_elliptic.H"
#include "xcontact_levelset_algorithm.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../solver_nonlin_nox/nox_nln_problem.H"
#include "../solver_nonlin_nox/nox_nln_globaldata.H"
#include "../solver_nonlin_nox/nox_nln_solver_factory.H"
#include <NOX_Solver_Generic.H>
#include <NOX_Abstract_Group.H>

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XCONTACT::LEVELSET::REINIT::Elliptic::Elliptic()
    : nox_params_(Teuchos::null),
      nox_problem_(Teuchos::null),
      grp_ptr_(Teuchos::null),
      ostatus_(Teuchos::null),
      jac_(Teuchos::null)
{
  /* left blank */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::REINIT::Elliptic::Setup()
{
  CheckInit();

  nox_params_ = Teuchos::rcp(new Teuchos::ParameterList());
  SetupNox(*nox_params_);

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::REINIT::Elliptic::PreSolve()
{
  if (Algorithm().Discretization()->Comm().MyPID() == 0)
  {
    std::cout << "*-------------------------------------------------------*\n";
    std::cout << "| SOLVE THE ELLIPTIC LEVEL-SET REINITIALIZATION PROBLEM |\n";
    std::cout << "*-------------------------------------------------------*\n";
    std::cout << std::flush;
  }

  Algorithm().SetReinitializationElementParameters();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::REINIT::Elliptic::Solve(const Epetra_Vector& phinp)
{
  Teuchos::RCP<NOX::Solver::Generic> nlnsolver = CreateNlnSolver(phinp);

  // solve the nonlinear problem
  nox_problem_->CheckFinalStatus(nlnsolver->solve());

  Group() = nlnsolver->getSolutionGroup();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Solver::Generic> XCONTACT::LEVELSET::REINIT::Elliptic::CreateNlnSolver(
    const Epetra_Vector& phinp)
{
  Algorithm().CreateActiveMaps(phinp);

  Teuchos::RCP<NOX::Epetra::Vector> soln = CreateSolutionVector(phinp);

  CreateJacobian(Algorithm().ActiveRowDofMap(), jac_);

  nox_problem_->Initialize(soln, jac_);

  Teuchos::RCP<NOX::Epetra::LinearSystem> linsystem = nox_problem_->CreateLinearSystem();

  grp_ptr_ = nox_problem_->CreateGroup(linsystem);

  return NOX::NLN::Solver::BuildSolver(
      grp_ptr_, ostatus_, Teuchos::null, nox_problem_->NlnGlobalDataPtr());
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::REINIT::Elliptic::CreateJacobian(
    const Epetra_Map& active_dofs, Teuchos::RCP<LINALG::SparseOperator>& jac)
{
  if (jac.is_null() or not jac->DomainMap().SameAs(active_dofs))
    jac = Teuchos::rcp(new LINALG::SparseMatrix(active_dofs, 9, false, true));
  else
    jac->Zero();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::Vector> XCONTACT::LEVELSET::REINIT::Elliptic::CreateSolutionVector(
    const Epetra_Vector& phinp_full) const
{
  Teuchos::RCP<Epetra_Vector> soln =
      Algorithm().ActiveDofMapExtractor().ExtractVector(phinp_full, XCONTACT::LEVELSET::active);

  // wrap and return
  return Teuchos::rcp(new NOX::Epetra::Vector(soln, NOX::Epetra::Vector::CreateView));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::REINIT::Elliptic::SetupNox(Teuchos::ParameterList& p_nox)
{
  /* Set necessary parameters */
  SetNlnSolverParams(p_nox);
  SetNlnSolverStatusTestParams(p_nox.sublist("Status Test"));
  SetNlnSolverPrintParams(p_nox.sublist("Printing"));

  // set linear solver
  std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>> linsolver;
  SetLinearSolver(linsolver, Algorithm().Solver());

  /* Set NOX::Epetra::Interface::Required */
  const Teuchos::RCP<NOX::Epetra::Interface::Required> ireq = Teuchos::rcp(this, false);

  /* Set NOX::Epetra::Interface::Jacobian */
  const Teuchos::RCP<NOX::Epetra::Interface::Jacobian> ijac = Teuchos::rcp(this, false);

  const Epetra_Comm& comm = Algorithm().Discretization()->Comm();
  Teuchos::RCP<NOX::NLN::GlobalData> nlnglobaldata =
      Teuchos::rcp(new NOX::NLN::GlobalData(comm, p_nox, linsolver, ireq, ijac));

  nox_problem_ = Teuchos::rcp(new NOX::NLN::Problem(nlnglobaldata));

  nox_problem_->CreateOuterStatusTest(ostatus_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::REINIT::Elliptic::SetNlnSolverParams(Teuchos::ParameterList& p_nox) const
{
  // nonlinear solver type
  p_nox.set("Nonlinear Solver", "Line Search Based");

  // solver options
  Teuchos::ParameterList& p_sol_opt = p_nox.sublist("Solver Options", false);
  p_sol_opt.set<std::string>("Status Test Check Type", "Complete");

  // direction method
  Teuchos::ParameterList& p_dir = p_nox.sublist("Direction", false);
  p_dir.set("Method", "Newton");
  Teuchos::ParameterList& p_newton = p_dir.sublist("Newton", false);
  Teuchos::ParameterList& p_lsolver = p_newton.sublist("Linear Solver", false);
  p_lsolver.set<bool>("Adaptive Control", false);

  // line search
  Teuchos::ParameterList& p_linesearch = p_nox.sublist("Line Search", false);
  p_linesearch.set("Method", "Full Step");

  // full step method
  Teuchos::ParameterList& p_fullstep = p_linesearch.sublist("Full Step", false);
  p_fullstep.set<double>("Full Step", 1.0);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::REINIT::Elliptic::SetNlnSolverStatusTestParams(
    Teuchos::ParameterList& p_status) const
{
  // get maximum number of iterations
  const int& max_iter = ReinitParams().get<int>("NUMSTEPSREINIT");

  // get tolerance
  const double& tol = ReinitParams().get<double>("CONVTOL_REINIT");

  Teuchos::ParameterList& p_ostatus = p_status.sublist("Outer Status Test", false);
  p_ostatus.set<std::string>("Test Type", "Combo");
  p_ostatus.set<std::string>("Combo Type", "OR");

  // === Stop if maximum number of iterations is reached ======================
  Teuchos::ParameterList& p_ostatus_0 = p_ostatus.sublist("Test 0", false);
  p_ostatus_0.set<std::string>("Test Type", "MaxIters");
  p_ostatus_0.set<int>("Maximum Iterations", max_iter);

  // === Stop if residual AND increment are small enough ======================
  Teuchos::ParameterList& p_ostatus_1 = p_ostatus.sublist("Test 1", false);
  p_ostatus_1.set<std::string>("Test Type", "Combo");
  p_ostatus_1.set<std::string>("Combo Type", "AND");

  // --- NormF test ( test the residual )
  Teuchos::ParameterList& p_ostatus_1_0 = p_ostatus_1.sublist("Test 0", false);
  p_ostatus_1_0.set<std::string>("Test Type", "NormF");
  p_ostatus_1_0.set<std::string>("Quantity Type", "LevelSet-Reinit");
  p_ostatus_1_0.set<std::string>("Tolerance Type", "Absolute");
  p_ostatus_1_0.set<double>("Tolerance", tol);
  p_ostatus_1_0.set<std::string>("Norm Type", "Two Norm");
  p_ostatus_1_0.set<std::string>("Scale Type", "Unscaled");

  // --- NormWRMS test ( test the increment )
  Teuchos::ParameterList& p_ostatus_1_1 = p_ostatus_1.sublist("Test 1", false);
  p_ostatus_1_1.set<double>("Alpha", 1.0);
  p_ostatus_1_1.set<double>("Beta", 0.5);
  p_ostatus_1_1.set<std::string>("Test Type", "NormWRMS");
  p_ostatus_1_1.set<std::string>("Quantity Type", "LevelSet-Reinit");
  p_ostatus_1_1.set<double>("Absolute Tolerance", tol);
  p_ostatus_1_1.set<double>("Relative Tolerance", 100 * tol);
  p_ostatus_1_1.set<double>("BDF Multiplier", 1.0);
  p_ostatus_1_1.set<std::string>("Disable Implicit Weighting", "Yes");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::REINIT::Elliptic::SetNlnSolverPrintParams(
    Teuchos::ParameterList& p_print) const
{
  p_print.set<bool>("Outer Iteration StatusTest", false);
  p_print.set<bool>("Inner Iteration", false);
  p_print.set<bool>("Outer Iteration", true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::LEVELSET::REINIT::Elliptic::SetLinearSolver(
    std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& linsolver,
    const Teuchos::RCP<LINALG::Solver>& lin_solver) const
{
  linsolver[NOX::NLN::sol_scatra] = lin_solver;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::ParameterList& XCONTACT::LEVELSET::REINIT::Elliptic::ReinitParams()
{
  return Algorithm().levelsetparams_->sublist("REINITIALIZATION", false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::ParameterList& XCONTACT::LEVELSET::REINIT::Elliptic::ReinitParams() const
{
  return Algorithm().levelsetparams_->sublist("REINITIALIZATION", false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::ParameterList& XCONTACT::LEVELSET::REINIT::Elliptic::NoxParams() { return *nox_params_; }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Teuchos::ParameterList& XCONTACT::LEVELSET::REINIT::Elliptic::NoxParams() const
{
  return *nox_params_;
}
