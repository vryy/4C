/*----------------------------------------------------------------------*/
/*! \file

\brief Base class for all FSI block preconditioning matrices

\level 1

*/
/*----------------------------------------------------------------------*/

#include <Epetra_Time.h>

#include "fsi_overlapprec.H"
#include "fsi_debugwriter.H"

#include "../drt_adapter/ad_ale_fsi.H"
#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "../drt_lib/drt_discret.H"

#include "../linalg/linalg_precond.H"
#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::BlockPreconditioningMatrix::BlockPreconditioningMatrix(
    Teuchos::RCP<UTILS::MonolithicDebugWriter> pcdbg, const LINALG::MultiMapExtractor& maps,
    ADAPTER::FSIStructureWrapper& structure, ADAPTER::Fluid& fluid, ADAPTER::AleFsiWrapper& ale,
    int symmetric, double omega, int iterations, double somega, int siterations, double fomega,
    int fiterations, double aomega, int aiterations, FILE* err)
    : LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(maps, maps, 81, false, true),
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
  fluidsolver_ = Teuchos::rcp(new LINALG::Preconditioner(fluid.LinearSolver()));

#ifndef BLOCKMATRIXMERGE
  structuresolver_ = Teuchos::rcp(new LINALG::Preconditioner(structure.LinearSolver()));

  constalesolver_ = ale.ConstPreconditioner();
  if (constalesolver_ == Teuchos::null)
    alesolver_ = Teuchos::rcp(new LINALG::Preconditioner(ale.LinearSolver()));
  else
    alesolver_ = constalesolver_;
#endif

  // check and fix ml nullspace if neccessary
  {
    LINALG::Solver& solver = *(structure.LinearSolver());
    const Epetra_Map& oldmap = *(structure.Discretization()->DofRowMap());
    const Epetra_Map& newmap = Matrix(0, 0).EpetraMatrix()->RowMap();
    solver.FixMLNullspace("Structure", oldmap, newmap, solver.Params());
  }
  {
    LINALG::Solver& solver = *(fluid.LinearSolver());
    const Epetra_Map& oldmap = *(fluid.DofRowMap());
    const Epetra_Map& newmap = Matrix(1, 1).EpetraMatrix()->RowMap();
    solver.FixMLNullspace("Fluid", oldmap, newmap, solver.Params());
  }
  {
    LINALG::Solver& solver = *(ale.LinearSolver());
    const Epetra_Map& oldmap = *(ale.Discretization()->DofRowMap());
    const Epetra_Map& newmap = Matrix(2, 2).EpetraMatrix()->RowMap();
    solver.FixMLNullspace("Ale", oldmap, newmap, solver.Params());
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
void FSI::BlockPreconditioningMatrix::MergeSolve(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
#ifdef BLOCKMATRIXMERGE
  const Epetra_Vector& x = Teuchos::dyn_cast<const Epetra_Vector>(X);
  Epetra_Vector& y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  fluidsolver_->Solve(
      sparse_->EpetraMatrix(), Teuchos::rcp(&y, false), Teuchos::rcp(new Epetra_Vector(x)), true);
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::BlockPreconditioningMatrix::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE

  // this is really evil :)
  sparse_ = Merge();
  fluidsolver_->Setup(sparse_->EpetraMatrix());

#else
  const LINALG::SparseMatrix& structInnerOp = Matrix(0, 0);
  const LINALG::SparseMatrix& fluidInnerOp = Matrix(1, 1);
  const LINALG::SparseMatrix& aleInnerOp = Matrix(2, 2);

  structuresolver_->Setup(structInnerOp.EpetraMatrix());
  fluidsolver_->Setup(fluidInnerOp.EpetraMatrix());
  if (constalesolver_ == Teuchos::null) alesolver_->Setup(aleInnerOp.EpetraMatrix());
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::BlockPreconditioningMatrix::LocalBlockRichardson(
    Teuchos::RCP<LINALG::Preconditioner> solver, const LINALG::SparseMatrix& innerOp,
    Teuchos::RCP<Epetra_Vector> x, Teuchos::RCP<Epetra_Vector> y, Teuchos::RCP<Epetra_Vector> tmpx,
    int iterations, double omega, FILE* err, const Epetra_Comm& comm)
{
  if (iterations > 0)
  {
    y->Scale(omega);
    Teuchos::RCP<Epetra_Vector> tmpy = Teuchos::rcp(new Epetra_Vector(y->Map()));
    if (err != NULL)
      if (comm.MyPID() == 0) fprintf(err, "    fluid richardson (%d,%f):", iterations, omega);
    for (int i = 0; i < iterations; ++i)
    {
      innerOp.EpetraMatrix()->Multiply(false, *y, *tmpx);
      tmpx->Update(1.0, *x, -1.0);

      if (err != NULL)
      {
        double n;
        tmpx->Norm2(&n);
        if (comm.MyPID() == 0) fprintf(err, " %e", n);
      }

      solver->Solve(innerOp.EpetraMatrix(), tmpy, tmpx, false);
      y->Update(omega, *tmpy, 1.0);
    }
    if (err != NULL)
      if (comm.MyPID() == 0) fprintf(err, "\n");
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::OverlappingBlockMatrix::OverlappingBlockMatrix(
    Teuchos::RCP<UTILS::MonolithicDebugWriter> pcdbg, const LINALG::MultiMapExtractor& maps,
    ADAPTER::FSIStructureWrapper& structure, ADAPTER::Fluid& fluid, ADAPTER::AleFsiWrapper& ale,
    bool structuresplit, int symmetric, double omega, int iterations, double somega,
    int siterations, double fomega, int fiterations, double aomega, int aiterations, FILE* err)
    : BlockPreconditioningMatrix(pcdbg, maps, structure, fluid, ale, symmetric, omega, iterations,
          somega, siterations, fomega, fiterations, aomega, aiterations, err),
      structuresplit_(structuresplit),
      structure_(structure),
      fluid_(fluid),
      ale_(ale)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE

  // this is really evil :)
  sparse_ = Merge();
  fluidsolver_->Setup(sparse_->EpetraMatrix());

#else
  const LINALG::SparseMatrix& structInnerOp = Matrix(0, 0);
  const LINALG::SparseMatrix& fluidInnerOp = Matrix(1, 1);
  const LINALG::SparseMatrix& aleInnerOp = Matrix(2, 2);

  Teuchos::RCP<LINALG::MapExtractor> fsidofmapex = Teuchos::null;
  Teuchos::RCP<Epetra_Map> irownodes = Teuchos::null;

  structuresolver_->Setup(structInnerOp.EpetraMatrix());
  fluidsolver_->Setup(fluidInnerOp.EpetraMatrix(), fsidofmapex, fluid_.Discretization(), irownodes,
      structuresplit_);
  if (constalesolver_ == Teuchos::null) alesolver_->Setup(aleInnerOp.EpetraMatrix());
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrix::SGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  const LINALG::SparseMatrix& structInnerOp = Matrix(0, 0);
  const LINALG::SparseMatrix& fluidInnerOp = Matrix(1, 1);
  const LINALG::SparseMatrix& fluidMeshOp = Matrix(1, 2);
  const LINALG::SparseMatrix& fluidBoundOp = Matrix(1, 0);
  const LINALG::SparseMatrix& aleInnerOp = Matrix(2, 2);

  // Extract vector blocks
  // RHS

  const Epetra_Vector& x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // initial guess

  Epetra_Vector& y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy = RangeExtractor().ExtractVector(y, 0);
  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().ExtractVector(y, 1);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().ExtractVector(y, 2);

  Teuchos::RCP<Epetra_Vector> sz = Teuchos::rcp(new Epetra_Vector(sy->Map()));
  Teuchos::RCP<Epetra_Vector> fz = Teuchos::rcp(new Epetra_Vector(fy->Map()));
  Teuchos::RCP<Epetra_Vector> az = Teuchos::rcp(new Epetra_Vector(ay->Map()));

  Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
  Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));

  // outer Richardson loop
  for (int run = 0; run < iterations_; ++run)
  {
    Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().ExtractVector(x, 0);
    Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x, 1);
    Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x, 2);

    // ----------------------------------------------------------------
    // lower GS

    {
      if (run > 0)
      {
        const LINALG::SparseMatrix& structBoundOp = Matrix(0, 1);

        structInnerOp.Multiply(false, *sy, *tmpsx);
        sx->Update(-1.0, *tmpsx, 1.0);
        structBoundOp.Multiply(false, *fy, *tmpsx);
        sx->Update(-1.0, *tmpsx, 1.0);
      }

      // Solve structure equations for sy with the rhs sx
      structuresolver_->Solve(structInnerOp.EpetraMatrix(), sz, sx, true);
      // do Richardson iteration
      LocalBlockRichardson(
          structuresolver_, structInnerOp, sx, sz, tmpsx, siterations_, somega_, err_, Comm());

      if (run > 0)
      {
        sy->Update(omega_, *sz, 1.0);
      }
      else
      {
        sy->Update(omega_, *sz, 0.0);
      }
    }

    {
      // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy

      if (run > 0)
      {
        aleInnerOp.Multiply(false, *ay, *tmpax);
        ax->Update(-1.0, *tmpax, 1.0);
      }

      if (structuresplit_)
      {
        if (run > 0)
        {
          const LINALG::SparseMatrix& aleBoundOp = Matrix(2, 1);
          aleBoundOp.Multiply(false, *fy, *tmpax);
          ax->Update(-1.0, *tmpax, 1.0);
        }
      }
      else
      {
        const LINALG::SparseMatrix& aleBoundOp = Matrix(2, 0);
        aleBoundOp.Multiply(false, *sy, *tmpax);
        ax->Update(-1.0, *tmpax, 1.0);
      }

      alesolver_->Solve(aleInnerOp.EpetraMatrix(), az, ax, true);
      // do Richardson iteration
      LocalBlockRichardson(
          alesolver_, aleInnerOp, ax, az, tmpax, aiterations_, aomega_, err_, Comm());

      if (run > 0)
      {
        ay->Update(omega_, *az, 1.0);
      }
      else
      {
        ay->Update(omega_, *az, 0.0);
      }
    }

    {
      // Solve fluid equations for fy with the rhs fx - F(I,Gamma) sy - F(Mesh) ay

      if (run > 0)
      {
        fluidInnerOp.Multiply(false, *fy, *tmpfx);
        fx->Update(-1.0, *tmpfx, 1.0);
      }

      fluidBoundOp.Multiply(false, *sy, *tmpfx);
      fx->Update(-1.0, *tmpfx, 1.0);
      fluidMeshOp.Multiply(false, *ay, *tmpfx);
      fx->Update(-1.0, *tmpfx, 1.0);
      fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(), fz, fx, true);

      LocalBlockRichardson(
          fluidsolver_, fluidInnerOp, fx, fz, tmpfx, fiterations_, fomega_, err_, Comm());

      if (run > 0)
      {
        fy->Update(omega_, *fz, 1.0);
      }
      else
      {
        fy->Update(omega_, *fz, 0.0);
      }
    }

    // ----------------------------------------------------------------
    // the symmetric part of the pc can be skipped

    if (symmetric_)
    {
      sx = DomainExtractor().ExtractVector(x, 0);
      fx = DomainExtractor().ExtractVector(x, 1);
      ax = DomainExtractor().ExtractVector(x, 2);

      // ----------------------------------------------------------------
      // upper GS

      {
        fluidInnerOp.Multiply(false, *fy, *tmpfx);
        fx->Update(-1.0, *tmpfx, 1.0);
        fluidBoundOp.Multiply(false, *sy, *tmpfx);
        fx->Update(-1.0, *tmpfx, 1.0);
        fluidMeshOp.Multiply(false, *ay, *tmpfx);
        fx->Update(-1.0, *tmpfx, 1.0);

        fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(), fz, fx, true);

        LocalBlockRichardson(
            fluidsolver_, fluidInnerOp, fx, fz, tmpfx, fiterations_, fomega_, err_, Comm());
        fy->Update(omega_, *fz, 1.0);
      }

      {
        aleInnerOp.Multiply(false, *ay, *tmpax);
        ax->Update(-1.0, *tmpax, 1.0);
        if (structuresplit_)
        {
          const LINALG::SparseMatrix& aleBoundOp = Matrix(2, 1);
          aleBoundOp.Multiply(false, *fy, *tmpax);
          ax->Update(-1.0, *tmpax, 1.0);
        }
        else
        {
          const LINALG::SparseMatrix& aleBoundOp = Matrix(2, 0);
          aleBoundOp.Multiply(false, *sy, *tmpax);
          ax->Update(-1.0, *tmpax, 1.0);
        }
        alesolver_->Solve(aleInnerOp.EpetraMatrix(), az, ax, true);

        LocalBlockRichardson(
            alesolver_, aleInnerOp, ax, az, tmpax, aiterations_, aomega_, err_, Comm());
        ay->Update(omega_, *az, 1.0);
      }

      {
        const LINALG::SparseMatrix& structBoundOp = Matrix(0, 1);

        structInnerOp.Multiply(false, *sy, *tmpsx);
        sx->Update(-1.0, *tmpsx, 1.0);
        structBoundOp.Multiply(false, *fy, *tmpsx);
        sx->Update(-1.0, *tmpsx, 1.0);

        structuresolver_->Solve(structInnerOp.EpetraMatrix(), sz, sx, true);

        LocalBlockRichardson(
            structuresolver_, structInnerOp, sx, sz, tmpsx, siterations_, somega_, err_, Comm());
        sy->Update(omega_, *sz, 1.0);
      }
    }
  }

  RangeExtractor().InsertVector(*sy, 0, y);
  RangeExtractor().InsertVector(*fy, 1, y);
  RangeExtractor().InsertVector(*ay, 2, y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::OverlappingBlockMatrix::Label() const { return "FSI::OverlappingBlockMatrix"; }
