/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of Baci's interface to Krylov solvers

\level 0

*/
/*---------------------------------------------------------------------*/

#include "4C_linear_solver_method_iterative.hpp"

#include "4C_linear_solver_amgnxn_preconditioner.hpp"
#include "4C_linear_solver_preconditioner_block.hpp"
#include "4C_linear_solver_preconditioner_ifpack.hpp"
#include "4C_linear_solver_preconditioner_krylovprojection.hpp"
#include "4C_linear_solver_preconditioner_ml.hpp"
#include "4C_linear_solver_preconditioner_muelu.hpp"
#include "4C_linear_solver_preconditioner_point.hpp"
#include "4C_utils_exceptions.hpp"

#include <BelosBiCGStabSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosLinearProblem.hpp>
#include <Epetra_Comm.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Map.h>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
CORE::LINEAR_SOLVER::IterativeSolver<MatrixType, VectorType>::IterativeSolver(
    const Epetra_Comm& comm, Teuchos::ParameterList& params)
    : comm_(comm), params_(params)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
void CORE::LINEAR_SOLVER::IterativeSolver<MatrixType, VectorType>::Setup(
    Teuchos::RCP<MatrixType> matrix, Teuchos::RCP<VectorType> x, Teuchos::RCP<VectorType> b,
    const bool refactor, const bool reset, Teuchos::RCP<CORE::LINALG::KrylovProjector> projector)
{
  // see whether operator is a Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(matrix);

  if (!Params().isSublist("Belos Parameters")) FOUR_C_THROW("Do not have belos parameter list");
  Teuchos::ParameterList& belist = Params().sublist("Belos Parameters");

  const int reuse = belist.get("reuse", 0);
  const bool create = !AllowReusePreconditioner(reuse, reset);
  if (create)
  {
    ncall_ = 0;
    preconditioner_ = CreatePreconditioner(belist, A != Teuchos::null, projector);
  }

  b_ = b;
  a_ = matrix;  // we cannot use A here, since it could be Teuchos::null (for blocked operators);
  x_ = x;

  preconditioner_->Setup(create, a_.get(), x_.get(), b_.get());
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
int CORE::LINEAR_SOLVER::IterativeSolver<MatrixType, VectorType>::Solve()
{
  Teuchos::ParameterList& belist = Params().sublist("Belos Parameters");

  auto problem = Teuchos::rcp(new Belos::LinearProblem<double, VectorType, MatrixType>(a_, x_, b_));

  if (preconditioner_ != Teuchos::null)
  {
    auto belosPrec = Teuchos::rcp(new Belos::EpetraPrecOp(preconditioner_->PrecOperator()));
    problem->setRightPrec(belosPrec);
  }

  const bool set = problem->setProblem();
  if (set == false)
    FOUR_C_THROW("CORE::LINEAR_SOLVER::BelosSolver: Iterative solver failed to set up correctly.");

  Teuchos::RCP<Belos::SolverManager<double, VectorType, MatrixType>> newSolver;
  std::string solverType = belist.get<std::string>("Solver Type");
  if (solverType == "GMRES")
    newSolver = Teuchos::rcp(new Belos::BlockGmresSolMgr<double, VectorType, MatrixType>(
        problem, Teuchos::rcp(&belist, false)));
  else if (solverType == "CG")
    newSolver = Teuchos::rcp(new Belos::BlockCGSolMgr<double, VectorType, MatrixType>(
        problem, Teuchos::rcp(&belist, false)));
  else if (solverType == "BiCGSTAB")
    newSolver = Teuchos::rcp(new Belos::BiCGStabSolMgr<double, VectorType, MatrixType>(
        problem, Teuchos::rcp(&belist, false)));
  else
    FOUR_C_THROW("CORE::LINEAR_SOLVER::BelosSolver: Unknown iterative solver solver type chosen.");

  Belos::ReturnType ret = newSolver->solve();

  numiters_ = newSolver->getNumIters();

  if (preconditioner_ != Teuchos::null) preconditioner_->Finish(a_.get(), x_.get(), b_.get());

  int my_error = 0;
  if (ret != Belos::Converged) my_error = 1;
  int glob_error = 0;
  comm_.SumAll(&my_error, &glob_error, 1);

  if (glob_error > 0 and this->comm_.MyPID() == 0)
    std::cout << std::endl
              << "CORE::LINEAR_SOLVER::BelosSolver: WARNING: Iterative solver did not converge!"
              << std::endl;

  ncall_ += 1;

  return 0;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
bool CORE::LINEAR_SOLVER::IterativeSolver<MatrixType, VectorType>::AllowReusePreconditioner(
    const int reuse, const bool reset)
{
  // first, check some parameters with information that has to be updated
  const Teuchos::ParameterList& linSysParams = Params().sublist("Belos Parameters");

  bool bAllowReuse = CheckReuseStatusOfActiveSet(linSysParams);

  const bool create = reset or not Ncall() or not reuse or (Ncall() % reuse) == 0;
  if (create) bAllowReuse = false;

  // here, each processor has its own local decision made
  // bAllowReuse = true -> preconditioner can be reused
  // bAllowReuse = false -> preconditioner has to be recomputed
  // If one or more processors decide that the preconditioner has to be recomputed
  // all of the processors have to recompute it

  // synchronize results of all processors
  // all processors have to do the same (either recompute preconditioner or allow reusing it)
  int nProc = comm_.NumProc();
  int lAllowReuse = bAllowReuse == true ? 1 : 0;
  int gAllowReuse = 0;
  comm_.SumAll(&lAllowReuse, &gAllowReuse, 1);

  return gAllowReuse == nProc;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
bool CORE::LINEAR_SOLVER::IterativeSolver<MatrixType, VectorType>::CheckReuseStatusOfActiveSet(
    const Teuchos::ParameterList& linSysParams)
{
  bool bAllowReuse = true;  // default: allow reuse of preconditioner

  if (linSysParams.isParameter("contact activeDofMap"))
  {
    auto ep_active_dof_map = linSysParams.get<Teuchos::RCP<Epetra_Map>>("contact activeDofMap");

    // Do we have history information available?
    if (active_dof_map_.is_null())
    {
      /* No history available.
       * This is the first application of the preconditioner. We cannot reuse it.
       */
      bAllowReuse = false;
    }
    else
    {
      /* History is available. We actually have to check for a change in the active set
       * by comparing the current map of active DOFs with the stored map of active DOFs
       * from the previous application of the preconditioner.
       */
      if (not ep_active_dof_map->PointSameAs(*active_dof_map_))
      {
        // Map of active nodes has changed -> force preconditioner to be rebuilt
        bAllowReuse = false;
      }
    }

    // Store current map of active slave DOFs for comparison in next preconditioner application
    active_dof_map_ = ep_active_dof_map;
  }

  return bAllowReuse;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
Teuchos::RCP<CORE::LINEAR_SOLVER::PreconditionerTypeBase>
CORE::LINEAR_SOLVER::IterativeSolver<MatrixType, VectorType>::CreatePreconditioner(
    Teuchos::ParameterList& solverlist, const bool isCrsMatrix,
    Teuchos::RCP<CORE::LINALG::KrylovProjector> projector)
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::Solver:  1.1)   CreatePreconditioner");

  Teuchos::RCP<CORE::LINEAR_SOLVER::PreconditionerTypeBase> preconditioner;

  if (isCrsMatrix)
  {
    // get type of preconditioner and build either Ifpack or ML
    // if we have an ifpack parameter list, we do ifpack
    // if we have an ml parameter list we do ml
    if (Params().isSublist("IFPACK Parameters"))
    {
      preconditioner = Teuchos::rcp(new CORE::LINEAR_SOLVER::IFPACKPreconditioner(
          Params().sublist("IFPACK Parameters"), solverlist));
    }
    else if (Params().isSublist("ML Parameters"))
    {
      preconditioner = Teuchos::rcp(
          new CORE::LINEAR_SOLVER::MLPreconditioner(Params().sublist("ML Parameters")));
    }
    else if (Params().isSublist("MueLu Parameters"))
    {
      preconditioner = Teuchos::rcp(
          new CORE::LINEAR_SOLVER::MueLuPreconditioner(Params().sublist("MueLu Parameters")));
    }
    else if (Params().isSublist("MueLu (BeamSolid) Parameters"))
    {
      preconditioner =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuBeamSolidBlockPreconditioner(Params()));
    }
    else
      FOUR_C_THROW(
          "CORE::LINEAR_SOLVER::IterativeSolver::CreatePreconditioner: Unknown preconditioner for "
          "iterative solver chosen.");

    // decide whether we do what kind of scaling
    std::string scaling = solverlist.get("scaling", "none");
    if (scaling == "none")
    {
    }
    else if (scaling == "infnorm")
    {
      preconditioner = Teuchos::rcp(new CORE::LINEAR_SOLVER::InfNormPreconditioner(preconditioner));
    }
    else if (scaling == "symmetric")
    {
      preconditioner = Teuchos::rcp(new CORE::LINEAR_SOLVER::SymDiagPreconditioner(preconditioner));
    }
    else
      FOUR_C_THROW(
          "CORE::LINEAR_SOLVER::IterativeSolver::CreatePreconditioner: Unknown type of scaling "
          "found in parameter list.");

    if (projector != Teuchos::null)
    {
      preconditioner = Teuchos::rcp(
          new CORE::LINEAR_SOLVER::KrylovProjectionPreconditioner(preconditioner, projector));
    }
  }
  else
  {
    // assume block matrix
    if (Params().isSublist("CheapSIMPLE Parameters"))
    {
      preconditioner = Teuchos::rcp(new CORE::LINEAR_SOLVER::SimplePreconditioner(Params()));
    }
    else if (Params().isSublist("BGS Parameters"))
    {
      preconditioner = Teuchos::rcp(
          new CORE::LINEAR_SOLVER::BGSPreconditioner(Params(), Params().sublist("BGS Parameters")));
    }
    else if (Params().isSublist("MueLu (Fluid) Parameters"))
    {
      preconditioner = Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuFluidBlockPreconditioner(
          Params().sublist("MueLu (Fluid) Parameters")));
    }
    else if (Params().isSublist("MueLu (TSI) Parameters"))
    {
      preconditioner = Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuTsiBlockPreconditioner(Params()));
    }
    else if (Params().isSublist("MueLu (Contact) Parameters"))
    {
      preconditioner =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuContactSpPreconditioner(Params()));
    }
    else if (Params().isSublist("MueLu (FSI) Parameters"))
    {
      preconditioner = Teuchos::rcp(new CORE::LINEAR_SOLVER::MueLuFsiBlockPreconditioner(Params()));
    }
    else if (Params().isSublist("AMGnxn Parameters"))
    {
      preconditioner = Teuchos::rcp(new CORE::LINEAR_SOLVER::AmGnxnPreconditioner(Params()));
    }
    else
      FOUR_C_THROW(
          "CORE::LINEAR_SOLVER::IterativeSolver::CreatePreconditioner: Unknown preconditioner for "
          "block iterative solver chosen.");
  }

  return preconditioner;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// explicit initialization
template class CORE::LINEAR_SOLVER::IterativeSolver<Epetra_Operator, Epetra_MultiVector>;

FOUR_C_NAMESPACE_CLOSE
