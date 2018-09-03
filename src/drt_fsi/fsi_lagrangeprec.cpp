

#include "fsi_lagrangeprec.H"
#include <Epetra_Time.h>

#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_ale_fsi.H"

#include "../drt_structure/stru_aux.H"

#include "../linalg/linalg_precond.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::LagrangianBlockMatrix::LagrangianBlockMatrix(const LINALG::MultiMapExtractor& maps,
    ADAPTER::FSIStructureWrapper& structure, ADAPTER::Fluid& fluid, ADAPTER::AleFsiWrapper& ale,
    int symmetric, double omega, int iterations, double somega, int siterations, double fomega,
    int fiterations, double aomega, int aiterations, FILE* err)
    : BlockPreconditioningMatrix(Teuchos::null, maps, structure, fluid, ale, symmetric, omega,
          iterations, somega, siterations, fomega, fiterations, aomega, aiterations, err),
      structure_(structure)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::LagrangianBlockMatrix::SGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  // Extract matrix blocks

  // const LINALG::SparseMatrix& S    = Matrix(0,0);
  const LINALG::SparseMatrix& F = Matrix(1, 1);
  const LINALG::SparseMatrix& Fg = Matrix(1, 2);
  const LINALG::SparseMatrix& Aii = Matrix(2, 2);
  const LINALG::SparseMatrix& Aig = Matrix(2, 1);

  const LINALG::SparseMatrix& CSF = Matrix(3, 0);
  const LINALG::SparseMatrix& CFS = Matrix(3, 1);

  const LINALG::SparseMatrix& CSFT = Matrix(0, 3);
  const LINALG::SparseMatrix& CFST = Matrix(1, 3);

  const LINALG::SparseMatrix& Sii = blocks_->Matrix(0, 0);
  const LINALG::SparseMatrix& Sig = blocks_->Matrix(0, 1);
  const LINALG::SparseMatrix& Sgi = blocks_->Matrix(1, 0);
  const LINALG::SparseMatrix& Sgg = blocks_->Matrix(1, 1);

  std::vector<Teuchos::RCP<const Epetra_Map>> innerstruct;
  innerstruct.push_back(structure_.Interface()->FSICondMap());
  innerstruct.push_back(structure_.Interface()->OtherMap());
  innerstruct.push_back(Teuchos::null);
  LINALG::MultiMapExtractor innerstructextract(FullRowMap(), innerstruct);

  // Extract vector blocks

  // RHS
  const Epetra_Vector& x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // initial guess
  Epetra_Vector& y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  // these are assumed to be zero
  Teuchos::RCP<Epetra_Vector> sy = RangeExtractor().ExtractVector(y, 0);
  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().ExtractVector(y, 1);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().ExtractVector(y, 2);
  Teuchos::RCP<Epetra_Vector> ly = RangeExtractor().ExtractVector(y, 3);

  Teuchos::RCP<Epetra_Vector> sgy =
      Teuchos::rcp(new Epetra_Vector(*structure_.Interface()->FSICondMap()));
  Teuchos::RCP<Epetra_Vector> siy =
      Teuchos::rcp(new Epetra_Vector(*structure_.Interface()->OtherMap()));

  Teuchos::RCP<Epetra_Vector> sz = Teuchos::rcp(new Epetra_Vector(sy->Map()));
  Teuchos::RCP<Epetra_Vector> fz = Teuchos::rcp(new Epetra_Vector(fy->Map()));
  Teuchos::RCP<Epetra_Vector> az = Teuchos::rcp(new Epetra_Vector(ay->Map()));

  Teuchos::RCP<Epetra_Vector> siz =
      Teuchos::rcp(new Epetra_Vector(*structure_.Interface()->OtherMap()));

  Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
  Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));
  Teuchos::RCP<Epetra_Vector> tmplx = Teuchos::rcp(new Epetra_Vector(DomainMap(3)));

  Teuchos::RCP<Epetra_Vector> tmpsgx = Teuchos::rcp(new Epetra_Vector(sgy->Map()));
  Teuchos::RCP<Epetra_Vector> tmpsix = Teuchos::rcp(new Epetra_Vector(siy->Map()));

  // block preconditioner

  if (symmetric_) dserror("symmetric Gauss-Seidel not implemented");

  // outer Richardson loop
  for (int run = 0; run < iterations_; ++run)
  {
    Teuchos::RCP<Epetra_Vector> sgx = innerstructextract.ExtractVector(x, 0);
    Teuchos::RCP<Epetra_Vector> six = innerstructextract.ExtractVector(x, 1);

    // Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().ExtractVector(x,0);
    Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x, 1);
    Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x, 2);
    Teuchos::RCP<Epetra_Vector> lx = DomainExtractor().ExtractVector(x, 3);

    // structure interface from fluid

    if (run > 0)
    {
      CFS.Multiply(false, *fy, *tmplx);
      lx->Update(-1.0, *tmplx, 1.0);

      // semi-solve. Assume CSF==I
      CSF.Multiply(true, *lx, *sy);
      structure_.Interface()->ExtractFSICondVector(sy, sgy);
    }

    // structure

    if (run > 0)
    {
      Sii.Multiply(false, *siy, *tmpsix);
      six->Update(-1.0, *tmpsix, 1.0);
      Sig.Multiply(false, *sgy, *tmpsix);
      six->Update(-1.0, *tmpsix, 1.0);
    }

    structuresolver_->Solve(Sii.EpetraMatrix(), siz, six, true);
    LocalBlockRichardson(
        structuresolver_, Sii, six, siz, tmpsix, siterations_, somega_, err_, Comm());

    if (run > 0)
    {
      siy->Update(omega_, *siz, 1.0);
    }
    else
    {
      siy->Update(omega_, *siz, 0.0);
    }

    // lagrange coupling

    if (run > 0)
    {
      Sgg.Multiply(false, *sgy, *tmpsgx);
      sgx->Update(-1.0, *tmpsgx, 1.0);
    }

    Sgi.Multiply(false, *siy, *tmpsgx);
    sgx->Update(-1.0, *tmpsgx, 1.0);

    // semi-solve. Assume CSFT==I
    structure_.Interface()->InsertFSICondVector(sgx, tmpsx);
    CSFT.Multiply(true, *tmpsx, *ly);

    // ale

    if (run > 0)
    {
      Aii.Multiply(false, *ay, *tmpax);
      ax->Update(-1.0, *tmpax, 1.0);
      Aig.Multiply(false, *fy, *tmpax);
      ax->Update(-1.0, *tmpax, 1.0);
    }

    alesolver_->Solve(Aii.EpetraMatrix(), az, ax, true);

    if (run > 0)
    {
      ay->Update(omega_, *az, 1.0);
    }
    else
    {
      ay->Update(omega_, *az, 0.0);
    }

    // fluid

    if (run > 0)
    {
      F.Multiply(false, *fy, *tmpfx);
      fx->Update(-1.0, *tmpfx, 1.0);
    }

    Fg.Multiply(false, *ay, *tmpfx);
    fx->Update(-1.0, *tmpfx, 1.0);
    CFST.Multiply(false, *ly, *tmpfx);
    fx->Update(-1.0, *tmpfx, 1.0);

    fluidsolver_->Solve(F.EpetraMatrix(), fz, fx, true);
    LocalBlockRichardson(fluidsolver_, F, fx, fz, tmpfx, fiterations_, fomega_, err_, Comm());

    if (run > 0)
    {
      fy->Update(omega_, *fz, 1.0);
    }
    else
    {
      fy->Update(omega_, *fz, 0.0);
    }
  }

  // build solution vector

  structure_.Interface()->InsertFSICondVector(sgy, sy);
  structure_.Interface()->InsertOtherVector(siy, sy);

  RangeExtractor().InsertVector(*sy, 0, y);
  RangeExtractor().InsertVector(*fy, 1, y);
  RangeExtractor().InsertVector(*ay, 2, y);
  // RangeExtractor().InsertVector(*ly,3,y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::LagrangianBlockMatrix::Label() const { return "FSI::LagrangianBlockMatrix"; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::LagrangianBlockMatrix::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE
  BlockPreconditioningMatrix::SetupPreconditioner();
#else
  const LINALG::SparseMatrix& structInnerOp = Matrix(0, 0);
  const LINALG::SparseMatrix& fluidInnerOp = Matrix(1, 1);
  const LINALG::SparseMatrix& aleInnerOp = Matrix(2, 2);

  blocks_ = structInnerOp.Split<LINALG::DefaultBlockMatrixStrategy>(
      *structure_.Interface(), *structure_.Interface());
  blocks_->Complete();

  LINALG::SparseMatrix& sii = blocks_->Matrix(0, 0);
  //   LINALG::SparseMatrix& sig = blocks_->Matrix(0,1);
  //   LINALG::SparseMatrix& sgi = blocks_->Matrix(1,0);
  //   LINALG::SparseMatrix& sgg = blocks_->Matrix(1,1);

  structuresolver_->Setup(sii.EpetraMatrix());
  fluidsolver_->Setup(fluidInnerOp.EpetraMatrix());
  if (constalesolver_ == Teuchos::null) alesolver_->Setup(aleInnerOp.EpetraMatrix());
#endif
}
