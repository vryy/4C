/*----------------------------------------------------------------------*/
/*! \file

\brief Base class for all FSI block preconditioning matrices

\level 1

*/
/*----------------------------------------------------------------------*/

#include "baci_fsi_overlapprec.hpp"

#include "baci_adapter_ale_fsi.hpp"
#include "baci_adapter_fld_fluid.hpp"
#include "baci_adapter_str_fsiwrapper.hpp"
#include "baci_fsi_debugwriter.hpp"
#include "baci_lib_discret.hpp"
#include "baci_linear_solver_method_linalg.hpp"
#include "baci_linear_solver_method_parameters.hpp"
#include "baci_linear_solver_preconditioner_linalg.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::BlockPreconditioningMatrix::BlockPreconditioningMatrix(
    Teuchos::RCP<UTILS::MonolithicDebugWriter> pcdbg, const CORE::LINALG::MultiMapExtractor& maps,
    ADAPTER::FSIStructureWrapper& structure, ADAPTER::Fluid& fluid, ADAPTER::AleFsiWrapper& ale,
    int symmetric, double omega, int iterations, double somega, int siterations, double fomega,
    int fiterations, double aomega, int aiterations, FILE* err)
    : CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
          maps, maps, 81, false, true),
      symmetric_(symmetric),
      omega_(omega),
      iterations_(iterations),
      somega_(somega),
      siterations_(siterations),
      fomega_(fomega),
      fiterations_(fiterations),
      aomega_(aomega),
      aiterations_(aiterations),
      err_(err),
      pcdbg_(pcdbg)
{
  fluidsolver_ = Teuchos::rcp(new CORE::LINALG::Preconditioner(fluid.LinearSolver()));

#ifndef BLOCKMATRIXMERGE
  structuresolver_ = Teuchos::rcp(new CORE::LINALG::Preconditioner(structure.LinearSolver()));
  alesolver_ = Teuchos::rcp(new CORE::LINALG::Preconditioner(ale.LinearSolver()));
#endif

  // check and fix ml nullspace if neccessary
  {
    CORE::LINALG::Solver& solver = *(structure.LinearSolver());
    const Epetra_Map& oldmap = *(structure.Discretization()->DofRowMap());
    const Epetra_Map& newmap = Matrix(0, 0).EpetraMatrix()->RowMap();
    CORE::LINEAR_SOLVER::Parameters::FixNullSpace("Structure", oldmap, newmap, solver.Params());
  }
  {
    CORE::LINALG::Solver& solver = *(fluid.LinearSolver());
    const Epetra_Map& oldmap = *(fluid.DofRowMap());
    const Epetra_Map& newmap = Matrix(1, 1).EpetraMatrix()->RowMap();
    CORE::LINEAR_SOLVER::Parameters::FixNullSpace("Fluid", oldmap, newmap, solver.Params());
  }
  {
    CORE::LINALG::Solver& solver = *(ale.LinearSolver());
    const Epetra_Map& oldmap = *(ale.Discretization()->DofRowMap());
    const Epetra_Map& newmap = Matrix(2, 2).EpetraMatrix()->RowMap();
    CORE::LINEAR_SOLVER::Parameters::FixNullSpace("Ale", oldmap, newmap, solver.Params());
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int FSI::BlockPreconditioningMatrix::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (UseTranspose()) dserror("no transpose preconditioning");

  Teuchos::RCP<Epetra_Vector> r;

  if (pcdbg_ != Teuchos::null)
  {
    pcdbg_->NewIteration();

    // X and Y are the same at this point (if we have been called by aztec!)
    Epetra_Vector& y = Teuchos::dyn_cast<Epetra_Vector>(Y);
    pcdbg_->WriteVector("x", Teuchos::rcp(&y, false));

    r = Teuchos::rcp(new Epetra_Vector(y.Map()));
    Apply(X, *r);
  }

#ifdef BLOCKMATRIXMERGE
  MergeSolve(X, Y);
#else
  SGS(X, Y);
#endif

  if (pcdbg_ != Teuchos::null)
  {
    Epetra_Vector& y = Teuchos::dyn_cast<Epetra_Vector>(Y);
    pcdbg_->WriteVector("y", Teuchos::rcp(&y, false));
    r->Update(-1, y, 1);
    pcdbg_->WriteVector("r", r);
  }

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::BlockPreconditioningMatrix::LocalBlockRichardson(
    Teuchos::RCP<CORE::LINALG::Preconditioner> solver, const CORE::LINALG::SparseMatrix& innerOp,
    Teuchos::RCP<Epetra_Vector> x, Teuchos::RCP<Epetra_Vector> y, Teuchos::RCP<Epetra_Vector> tmpx,
    int iterations, double omega, FILE* err, const Epetra_Comm& comm)
{
  if (iterations > 0)
  {
    y->Scale(omega);
    Teuchos::RCP<Epetra_Vector> tmpy = Teuchos::rcp(new Epetra_Vector(y->Map()));
    if (err != nullptr)
      if (comm.MyPID() == 0) fprintf(err, "    fluid richardson (%d,%f):", iterations, omega);
    for (int i = 0; i < iterations; ++i)
    {
      innerOp.EpetraMatrix()->Multiply(false, *y, *tmpx);
      tmpx->Update(1.0, *x, -1.0);

      if (err != nullptr)
      {
        double n;
        tmpx->Norm2(&n);
        if (comm.MyPID() == 0) fprintf(err, " %e", n);
      }

      solver->Solve(innerOp.EpetraMatrix(), tmpy, tmpx, false);
      y->Update(omega, *tmpy, 1.0);
    }
    if (err != nullptr)
      if (comm.MyPID() == 0) fprintf(err, "\n");
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::OverlappingBlockMatrix::OverlappingBlockMatrix(
    Teuchos::RCP<UTILS::MonolithicDebugWriter> pcdbg, const CORE::LINALG::MultiMapExtractor& maps,
    ADAPTER::FSIStructureWrapper& structure, ADAPTER::Fluid& fluid, ADAPTER::AleFsiWrapper& ale,
    bool structuresplit, int symmetric, double omega, int iterations, double somega,
    int siterations, double fomega, int fiterations, double aomega, int aiterations)
    : BlockPreconditioningMatrix(pcdbg, maps, structure, fluid, ale, symmetric, omega, iterations,
          somega, siterations, fomega, fiterations, aomega, aiterations),
      structuresplit_(structuresplit),
      structure_(structure),
      fluid_(fluid),
      ale_(ale)
{
}

BACI_NAMESPACE_CLOSE
