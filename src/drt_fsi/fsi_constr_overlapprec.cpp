/*----------------------------------------------------------------------*/
/*!
\file fsi_constr_overlapprec.cpp

\level 1

\maintainer Matthias Mayr

\brief Preconditioner for FSI problems with additional constraints
*/
/*----------------------------------------------------------------------*/

#include <Epetra_Time.h>

#include "fsi_constr_overlapprec.H"
#include "fsi_debugwriter.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid.H"
#include "../linalg/linalg_precond.H"
#include "../linalg/linalg_solver.H"

// /*----------------------------------------------------------------------*
//  *----------------------------------------------------------------------*/
FSI::ConstrOverlappingBlockMatrix::ConstrOverlappingBlockMatrix(
    const LINALG::MultiMapExtractor& maps, ADAPTER::FSIStructureWrapper& structure,
    ADAPTER::Fluid& fluid, ADAPTER::AleFsiWrapper& ale, bool structuresplit, int symmetric,
    double omega, int iterations, double somega, int siterations, double fomega, int fiterations,
    double aomega, int aiterations, FILE* err)
    : FSI::OverlappingBlockMatrix(Teuchos::null, maps, structure, fluid, ale, structuresplit,
          symmetric, omega, iterations, somega, siterations, fomega, fiterations, aomega,
          aiterations, err)
{
  // determine map of all dofs not related to constraint

  std::vector<Teuchos::RCP<const Epetra_Map>> fsimaps;
  fsimaps.push_back(maps.Map(0));
  fsimaps.push_back(maps.Map(1));
  fsimaps.push_back(maps.Map(2));
  overallfsimap_ = LINALG::MultiMapExtractor::MergeMaps(fsimaps);
  fsiextractor_ = LINALG::MultiMapExtractor(*overallfsimap_, fsimaps);

  // stuff needed for SIMPLE preconditioner -> this needs to be read
  // in from the input file one day!
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  alpha_ = fsidyn.sublist("CONSTRAINT").get<double>("ALPHA");
  simpleiter_ = fsidyn.sublist("CONSTRAINT").get<int>("SIMPLEITER");
  prec_ = DRT::INPUT::IntegralValue<INPAR::FSI::PrecConstr>(
      fsidyn.sublist("CONSTRAINT"), "PRECONDITIONER");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::ConstrOverlappingBlockMatrix::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE

  // this is really evil
  // note that this only works with UMFPACK as fluid solver (saddle
  // point problem!)

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
void FSI::ConstrOverlappingBlockMatrix::SGS(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  const LINALG::SparseMatrix& StructInnerOp = Matrix(0, 0);
  const LINALG::SparseMatrix& StructBoundOp = Matrix(0, 1);
  const LINALG::SparseMatrix& StructConOp = Matrix(0, 3);
  const LINALG::SparseMatrix& FluidBoundOp = Matrix(1, 0);
  const LINALG::SparseMatrix& FluidInnerOp = Matrix(1, 1);
  const LINALG::SparseMatrix& FluidMeshOp = Matrix(1, 2);
  const LINALG::SparseMatrix& FluidConOp = Matrix(1, 3);
  const LINALG::SparseMatrix& AleInnerOp = Matrix(2, 2);
  const LINALG::SparseMatrix& ConStructOp = Matrix(3, 0);
  const LINALG::SparseMatrix& ConFluidOp = Matrix(3, 1);


  // Extract vector blocks
  // RHS

  const Epetra_Vector& x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // initial guess

  Epetra_Vector& y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy = RangeExtractor().ExtractVector(y, 0);
  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().ExtractVector(y, 1);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().ExtractVector(y, 2);
  Teuchos::RCP<Epetra_Vector> cy = RangeExtractor().ExtractVector(y, 3);

  Teuchos::RCP<Epetra_Vector> sz = Teuchos::rcp(new Epetra_Vector(sy->Map()));
  Teuchos::RCP<Epetra_Vector> fz = Teuchos::rcp(new Epetra_Vector(fy->Map()));
  Teuchos::RCP<Epetra_Vector> az = Teuchos::rcp(new Epetra_Vector(ay->Map()));

  Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));
  Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(2)));

  for (int outerrun = 0; outerrun < simpleiter_; ++outerrun)
  {
    // -------------------------------------------------------------------
    // intermediate fsi dofs: u_(n+1/2) = F^(-1)(f-B^T*lambda_n)
    // -------------------------------------------------------------------

    // inner Richardson loop (FSI block)

    for (int run = 0; run < iterations_; ++run)
    {
      Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().ExtractVector(x, 0);
      Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x, 1);
      Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x, 2);

      // ----------------------------------------------------------------
      // lower GS


      // Structure
      {
        if (run > 0 or outerrun > 0)
        {
          StructConOp.Multiply(false, *cy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);

          StructInnerOp.Multiply(false, *sy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);
          StructBoundOp.Multiply(false, *fy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);
        }

        // Solve structure equations for sy with the rhs sx
        structuresolver_->Solve(StructInnerOp.EpetraMatrix(), sz, sx, true);

        // do Richardson iteration
        LocalBlockRichardson(
            structuresolver_, StructInnerOp, sx, sz, tmpsx, siterations_, somega_, err_, Comm());

        if (run > 0 or outerrun > 0)
        {
          sy->Update(omega_, *sz, 1.0);
        }
        else
        {
          sy->Update(omega_, *sz, 0.0);
        }
      }

      // Ale
      {
        // Solve ale equations for ay with the rhs ax - A(I,Gamma) sy

        if (run > 0 or outerrun > 0)
        {
          AleInnerOp.Multiply(false, *ay, *tmpax);
          ax->Update(-1.0, *tmpax, 1.0);
        }
        if (structuresplit_)
        {
          if (run > 0 or outerrun > 0)
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

        alesolver_->Solve(AleInnerOp.EpetraMatrix(), az, ax, true);

        // do Richardson iteration
        LocalBlockRichardson(
            alesolver_, AleInnerOp, ax, az, tmpax, aiterations_, aomega_, err_, Comm());

        if (run > 0 or outerrun > 0)
        {
          ay->Update(omega_, *az, 1.0);
        }
        else
        {
          ay->Update(omega_, *az, 0.0);
        }
      }

      // Fluid
      {
        // Solve fluid equations for fy with the rhs fx - F(I,Gamma) sy - F(Mesh) ay - F(Constr) cy

        if (run > 0 or outerrun > 0)
        {
          FluidInnerOp.Multiply(false, *fy, *tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);

          FluidConOp.Multiply(false, *cy, *tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);
        }

        FluidBoundOp.Multiply(false, *sy, *tmpfx);
        fx->Update(-1.0, *tmpfx, 1.0);
        FluidMeshOp.Multiply(false, *ay, *tmpfx);
        fx->Update(-1.0, *tmpfx, 1.0);


        fluidsolver_->Solve(FluidInnerOp.EpetraMatrix(), fz, fx, true);

        LocalBlockRichardson(
            fluidsolver_, FluidInnerOp, fx, fz, tmpfx, fiterations_, fomega_, err_, Comm());

        if (run > 0 or outerrun > 0)
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
          FluidInnerOp.Multiply(false, *fy, *tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);
          FluidBoundOp.Multiply(false, *sy, *tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);
          FluidMeshOp.Multiply(false, *ay, *tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);
          FluidConOp.Multiply(false, *cy, *tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);

          fluidsolver_->Solve(FluidInnerOp.EpetraMatrix(), fz, fx, true);

          LocalBlockRichardson(
              fluidsolver_, FluidInnerOp, fx, fz, tmpfx, fiterations_, fomega_, err_, Comm());
          fy->Update(omega_, *fz, 1.0);
        }

        {
          AleInnerOp.Multiply(false, *ay, *tmpax);
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

          alesolver_->Solve(AleInnerOp.EpetraMatrix(), az, ax, true);

          LocalBlockRichardson(
              alesolver_, AleInnerOp, ax, az, tmpax, aiterations_, aomega_, err_, Comm());

          ay->Update(omega_, *az, 1.0);
        }

        {
          StructInnerOp.Multiply(false, *sy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);
          StructBoundOp.Multiply(false, *fy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);
          StructConOp.Multiply(false, *cy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);


          structuresolver_->Solve(StructInnerOp.EpetraMatrix(), sz, sx, true);

          LocalBlockRichardson(
              structuresolver_, StructInnerOp, sx, sz, tmpsx, siterations_, somega_, err_, Comm());
          sy->Update(omega_, *sz, 1.0);
        }
      }
    }

    // -----------------------------------------------------------------------
    // intermediate constraint dofs: Dlambda~ = S^(-1) * (cx - B^ * u_(n+1/2))
    // -----------------------------------------------------------------------

    // cx - B^ * u_(n+1/2)

    Teuchos::RCP<Epetra_Vector> cx = DomainExtractor().ExtractVector(x, 3);

    Epetra_Vector inter(cx->Map());
    ConStructOp.Multiply(false, *sy, inter);
    cx->Update(-1.0, inter, 1.0);
    ConFluidOp.Multiply(false, *fy, inter);
    cx->Update(-1.0, inter, 1.0);


    LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> invDiag(
        fsiextractor_, fsiextractor_, 1, false, true);

    if (prec_ == INPAR::FSI::Simple)
    {
      // D^{-1} = diag(A(0,0))^{-1}

      Epetra_Vector structDiagVec(StructInnerOp.RowMap(), false);
      StructInnerOp.ExtractDiagonalCopy(structDiagVec);
      int err = structDiagVec.Reciprocal(structDiagVec);
      if (err) dserror("Epetra_MultiVector::Reciprocal (structure matrix) returned %d", err);
      LINALG::SparseMatrix invstructDiag(structDiagVec);

      Epetra_Vector fluidDiagVec(FluidInnerOp.RowMap(), false);
      FluidInnerOp.ExtractDiagonalCopy(fluidDiagVec);
      err = fluidDiagVec.Reciprocal(fluidDiagVec);
      if (err) dserror("Epetra_MultiVector::Reciprocal (fluid matrix) returned %d", err);
      LINALG::SparseMatrix invfluidDiag(fluidDiagVec);

      Epetra_Vector aleDiagVec(AleInnerOp.RowMap(), false);
      AleInnerOp.ExtractDiagonalCopy(aleDiagVec);
      err = aleDiagVec.Reciprocal(aleDiagVec);
      if (err) dserror("Epetra_MultiVector::Reciprocal (ale matrix) returned %d", err);
      LINALG::SparseMatrix invaleDiag(aleDiagVec);

      invDiag.Assign(0, 0, LINALG::View, invstructDiag);
      invDiag.Assign(1, 1, LINALG::View, invfluidDiag);
      invDiag.Assign(2, 2, LINALG::View, invaleDiag);
    }
    else if (prec_ == INPAR::FSI::Simplec)
    {
      // D^{-1} = sum(abs(A(0,0)))^{-1}

      Epetra_Vector structDiagVec(StructInnerOp.RowMap(), false);
      int err = StructInnerOp.EpetraMatrix()->InvRowSums(structDiagVec);
      if (err) dserror("Epetra_CrsMatrix::InvRowSums (structure matrix) returned %d", err);
      LINALG::SparseMatrix invstructDiag(structDiagVec);

      Epetra_Vector fluidDiagVec(FluidInnerOp.RowMap(), false);
      err = FluidInnerOp.EpetraMatrix()->InvRowSums(fluidDiagVec);
      if (err) dserror("Epetra_CrsMatrix::InvRowSums (fluid matrix) returned %d", err);
      LINALG::SparseMatrix invfluidDiag(fluidDiagVec);

      Epetra_Vector aleDiagVec(AleInnerOp.RowMap(), false);
      err = AleInnerOp.EpetraMatrix()->InvRowSums(aleDiagVec);
      if (err) dserror("Epetra_CrsMatrix::InvRowSums (ale matrix) returned %d", err);
      LINALG::SparseMatrix invaleDiag(aleDiagVec);

      invDiag.Assign(0, 0, LINALG::View, invstructDiag);
      invDiag.Assign(1, 1, LINALG::View, invfluidDiag);
      invDiag.Assign(2, 2, LINALG::View, invaleDiag);
    }
    else
      dserror("Unknown type of preconditioner for constraint fsi system");

    invDiag.Complete();


    // S = - B^ * D^{-1} * B^T

    Teuchos::RCP<LINALG::SparseMatrix> temps =
        LINALG::Multiply(ConStructOp, false, invDiag.Matrix(0, 0), false);
    Teuchos::RCP<LINALG::SparseMatrix> interconA =
        LINALG::Multiply(*temps, false, StructConOp, false, false);

    temps = LINALG::Multiply(ConFluidOp, false, invDiag.Matrix(1, 1), false);
    Teuchos::RCP<LINALG::SparseMatrix> tempss = LINALG::Multiply(*temps, false, FluidConOp, false);
    interconA->Add(*tempss, false, 1.0, 1.0);
    interconA->Complete(StructConOp.DomainMap(), ConStructOp.RangeMap());
    interconA->Scale(-1.0);

    Teuchos::RCP<Teuchos::ParameterList> constrsolvparams =
        Teuchos::rcp(new Teuchos::ParameterList);
    constrsolvparams->set("solver", "umfpack");
    Teuchos::RCP<Epetra_Vector> interconsol =
        Teuchos::rcp(new Epetra_Vector(ConStructOp.RangeMap()));
    Teuchos::RCP<LINALG::Solver> ConstraintSolver = Teuchos::rcp(new LINALG::Solver(
        constrsolvparams, interconA->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
    ConstraintSolver->Solve(interconA->EpetraOperator(), interconsol, cx, true, true);

    // -------------------------------------------------------------------
    // update of all dofs
    // -------------------------------------------------------------------

    if (outerrun > 0)
      cy->Update(alpha_, *interconsol, 1.0);
    else
      cy->Update(alpha_, *interconsol, 0.0);

    Teuchos::RCP<Epetra_Vector> temp1;
    Teuchos::RCP<Epetra_Vector> temp2;

    temp1 = Teuchos::rcp(new Epetra_Vector(sy->Map()));
    temp2 = Teuchos::rcp(new Epetra_Vector(sy->Map()));
    StructConOp.Multiply(false, *interconsol, *temp1);
    temp1->Scale(alpha_);
    invDiag.Matrix(0, 0).Multiply(false, *temp1, *temp2);
    sy->Update(-1.0, *temp2, 1.0);

    temp1 = Teuchos::rcp(new Epetra_Vector(fy->Map()));
    temp2 = Teuchos::rcp(new Epetra_Vector(fy->Map()));
    FluidConOp.Multiply(false, *interconsol, *temp1);
    temp1->Scale(alpha_);
    invDiag.Matrix(1, 1).Multiply(false, *temp1, *temp2);
    fy->Update(-1.0, *temp2, 1.0);
  }

  RangeExtractor().InsertVector(*sy, 0, y);
  RangeExtractor().InsertVector(*fy, 1, y);
  RangeExtractor().InsertVector(*ay, 2, y);
  RangeExtractor().InsertVector(*cy, 3, y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::ConstrOverlappingBlockMatrix::Label() const
{
  return "FSI::ConstrOverlappingBlockMatrix";
}
