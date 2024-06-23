/*----------------------------------------------------------------------*/
/*! \file
\brief BGS preconditioner for volume-coupled FSI


\level 2
*/
/*----------------------------------------------------------------------*/


#include "4C_fsi_lung_overlapprec.hpp"

#include "4C_adapter_fld_fluid.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_multiply.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_preconditioner_linalg.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

// /*----------------------------------------------------------------------*
//  *----------------------------------------------------------------------*/
FSI::LungOverlappingBlockMatrix::LungOverlappingBlockMatrix(
    const Core::LinAlg::MultiMapExtractor& maps, Adapter::FSIStructureWrapper& structure,
    Adapter::Fluid& fluid, Adapter::AleFsiWrapper& ale, bool structuresplit, int symmetric,
    double omega, int iterations, double somega, int siterations, double fomega, int fiterations,
    double aomega, int aiterations)
    : OverlappingBlockMatrix(Teuchos::null, maps, structure, fluid, ale, structuresplit, symmetric,
          omega, iterations, somega, siterations, fomega, fiterations, aomega, aiterations)
{
  // determine map of all dofs not related to constraint

  std::vector<Teuchos::RCP<const Epetra_Map>> fsimaps;
  fsimaps.push_back(maps.Map(0));
  fsimaps.push_back(maps.Map(1));
  fsimaps.push_back(maps.Map(2));
  overallfsimap_ = Core::LinAlg::MultiMapExtractor::merge_maps(fsimaps);
  fsiextractor_ = Core::LinAlg::MultiMapExtractor(*overallfsimap_, fsimaps);

  StructSchur_ = Teuchos::rcp(new LungSchurComplement());
  FluidSchur_ = Teuchos::rcp(new LungSchurComplement());

  // stuff needed for SIMPLE preconditioner
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();
  alpha_ = fsidyn.sublist("CONSTRAINT").get<double>("ALPHA");
  simpleiter_ = fsidyn.sublist("CONSTRAINT").get<int>("SIMPLEITER");
  prec_ = Core::UTILS::IntegralValue<Inpar::FSI::PrecConstr>(
      fsidyn.sublist("CONSTRAINT"), "PRECONDITIONER");

  Teuchos::ParameterList constrsolvparams;
  Core::UTILS::AddEnumClassToParameterList<Core::LinearSolver::SolverType>(
      "SOLVER", Core::LinearSolver::SolverType::umfpack, constrsolvparams);
  constraintsolver_ = Teuchos::rcp(new Core::LinAlg::Solver(constrsolvparams, maps.Map(0)->Comm(),
      Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY")));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::LungOverlappingBlockMatrix::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE

  // this is really evil :)
  // note that this only works with UMFPACK as fluid solver (saddle
  // point problem!)

  sparse_ = Merge();
  fluidsolver_->setup(sparse_->EpetraMatrix());

#else

  const Core::LinAlg::SparseMatrix& structInnerOp = Matrix(0, 0);
  const Core::LinAlg::SparseMatrix& fluidInnerOp = Matrix(1, 1);
  const Core::LinAlg::SparseMatrix& aleInnerOp = Matrix(2, 2);

  Teuchos::RCP<Core::LinAlg::MapExtractor> fsidofmapex = Teuchos::null;
  Teuchos::RCP<Epetra_Map> irownodes = Teuchos::null;

  structuresolver_->setup(structInnerOp.EpetraMatrix());
  fluidsolver_->setup(fluidInnerOp.EpetraMatrix(), fsidofmapex, fluid_.discretization(), irownodes,
      structuresplit_);
  if (constalesolver_ == Teuchos::null) alesolver_->setup(aleInnerOp.EpetraMatrix());


  // We can compute the schur complement only once
  invDiag_ =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          fsiextractor_, fsiextractor_, 1, false, true));

  const Core::LinAlg::SparseMatrix& StructInnerOp = Matrix(0, 0);
  const Core::LinAlg::SparseMatrix& StructConOp = Matrix(0, 3);
  const Core::LinAlg::SparseMatrix& FluidInnerOp = Matrix(1, 1);
  const Core::LinAlg::SparseMatrix& AleInnerOp = Matrix(2, 2);
  const Core::LinAlg::SparseMatrix& ConStructOp = Matrix(3, 0);

  if (prec_ == Inpar::FSI::Simple)
  {
    // D^{-1} = diag(A(0,0))^{-1}

    Epetra_Vector structDiagVec(StructInnerOp.RowMap(), false);
    StructInnerOp.ExtractDiagonalCopy(structDiagVec);
    int err = structDiagVec.Reciprocal(structDiagVec);
    if (err) FOUR_C_THROW("Epetra_MultiVector::Reciprocal (structure matrix) returned %d", err);
    Core::LinAlg::SparseMatrix invstructDiag(structDiagVec);

    Epetra_Vector fluidDiagVec(FluidInnerOp.RowMap(), false);
    FluidInnerOp.ExtractDiagonalCopy(fluidDiagVec);
    err = fluidDiagVec.Reciprocal(fluidDiagVec);
    if (err) FOUR_C_THROW("Epetra_MultiVector::Reciprocal (fluid matrix) returned %d", err);
    Core::LinAlg::SparseMatrix invfluidDiag(fluidDiagVec);

    Epetra_Vector aleDiagVec(AleInnerOp.RowMap(), false);
    AleInnerOp.ExtractDiagonalCopy(aleDiagVec);
    err = aleDiagVec.Reciprocal(aleDiagVec);
    if (err) FOUR_C_THROW("Epetra_MultiVector::Reciprocal (ale matrix) returned %d", err);
    Core::LinAlg::SparseMatrix invaleDiag(aleDiagVec);

    invDiag_->Assign(0, 0, Core::LinAlg::View, invstructDiag);
    invDiag_->Assign(1, 1, Core::LinAlg::View, invfluidDiag);
    invDiag_->Assign(2, 2, Core::LinAlg::View, invaleDiag);
  }
  else if (prec_ == Inpar::FSI::Simplec)
  {
    // D^{-1} = sum(abs(A(0,0)))^{-1}

    Epetra_Vector structDiagVec(StructInnerOp.RowMap(), false);
    int err = StructInnerOp.EpetraMatrix()->InvRowSums(structDiagVec);
    if (err) FOUR_C_THROW("Epetra_CrsMatrix::InvRowSums (structure matrix) returned %d", err);
    Core::LinAlg::SparseMatrix invstructDiag(structDiagVec);

    Epetra_Vector fluidDiagVec(FluidInnerOp.RowMap(), false);
    err = FluidInnerOp.EpetraMatrix()->InvRowSums(fluidDiagVec);
    if (err) FOUR_C_THROW("Epetra_CrsMatrix::InvRowSums (fluid matrix) returned %d", err);
    Core::LinAlg::SparseMatrix invfluidDiag(fluidDiagVec);

    Epetra_Vector aleDiagVec(AleInnerOp.RowMap(), false);
    err = AleInnerOp.EpetraMatrix()->InvRowSums(aleDiagVec);
    if (err) FOUR_C_THROW("Epetra_CrsMatrix::InvRowSums (ale matrix) returned %d", err);
    Core::LinAlg::SparseMatrix invaleDiag(aleDiagVec);

    invDiag_->Assign(0, 0, Core::LinAlg::View, invstructDiag);
    invDiag_->Assign(1, 1, Core::LinAlg::View, invfluidDiag);
    invDiag_->Assign(2, 2, Core::LinAlg::View, invaleDiag);
  }
  else
    FOUR_C_THROW("Unknown type of preconditioner for constraint fsi system");

  invDiag_->Complete();


  // S = - B^ * D^{-1} * B^T

  interconA_ = StructSchur_->CalculateSchur(Matrix(3, 0), invDiag_->Matrix(0, 0), Matrix(0, 3));
  Teuchos::RCP<Core::LinAlg::SparseMatrix> temp =
      FluidSchur_->CalculateSchur(Matrix(3, 1), invDiag_->Matrix(1, 1), Matrix(1, 3));

  interconA_->Add(*temp, false, 1.0, 1.0);

  interconA_->Complete(StructConOp.DomainMap(), ConStructOp.RangeMap());
  interconA_->Scale(-1.0);



#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::LungOverlappingBlockMatrix::sgs(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  const Core::LinAlg::SparseMatrix& StructInnerOp = Matrix(0, 0);
  const Core::LinAlg::SparseMatrix& StructBoundOp = Matrix(0, 1);
  const Core::LinAlg::SparseMatrix& StructConOp = Matrix(0, 3);
  const Core::LinAlg::SparseMatrix& FluidBoundOp = Matrix(1, 0);
  const Core::LinAlg::SparseMatrix& FluidInnerOp = Matrix(1, 1);
  const Core::LinAlg::SparseMatrix& FluidMeshOp = Matrix(1, 2);
  const Core::LinAlg::SparseMatrix& FluidConOp = Matrix(1, 3);
  const Core::LinAlg::SparseMatrix& AleInnerOp = Matrix(2, 2);
  const Core::LinAlg::SparseMatrix& ConStructOp = Matrix(3, 0);
  const Core::LinAlg::SparseMatrix& ConFluidOp = Matrix(3, 1);


  // Extract vector blocks
  // RHS

  const Epetra_Vector& x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // initial guess

  Epetra_Vector& y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy = RangeExtractor().extract_vector(y, 0);
  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().extract_vector(y, 1);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().extract_vector(y, 2);
  Teuchos::RCP<Epetra_Vector> cy = RangeExtractor().extract_vector(y, 3);

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
      Teuchos::RCP<Epetra_Vector> sx = DomainExtractor().extract_vector(x, 0);
      Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().extract_vector(x, 1);
      Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().extract_vector(x, 2);

      // ----------------------------------------------------------------
      // lower GS


      // Structure
      {
        if (run > 0 or outerrun > 0)
        {
          StructInnerOp.Multiply(false, *sy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);
          StructBoundOp.Multiply(false, *fy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);
          StructConOp.Multiply(false, *cy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);
        }

        // Solve structure equations for sy with the rhs sx
        structuresolver_->Solve(StructInnerOp.EpetraMatrix(), sz, sx, true);

        // do Richardson iteration
        local_block_richardson(
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
            const Core::LinAlg::SparseMatrix& aleFBoundOp = Matrix(2, 1);
            aleFBoundOp.Multiply(false, *fy, *tmpax);
            ax->Update(-1.0, *tmpax, 1.0);
          }

          const Core::LinAlg::SparseMatrix& aleSBoundOp = Matrix(2, 0);
          aleSBoundOp.Multiply(false, *sy, *tmpax);
          ax->Update(-1.0, *tmpax, 1.0);
        }
        else
        {
          const Core::LinAlg::SparseMatrix& aleFBoundOp = Matrix(2, 0);
          aleFBoundOp.Multiply(false, *sy, *tmpax);
          ax->Update(-1.0, *tmpax, 1.0);
        }

        alesolver_->Solve(AleInnerOp.EpetraMatrix(), az, ax, true);

        // do Richardson iteration
        local_block_richardson(
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

        local_block_richardson(
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
    }

    // -----------------------------------------------------------------------
    // intermediate constraint dofs: Dlambda~ = S^(-1) * (cx - B^ * u_(n+1/2))
    // -----------------------------------------------------------------------

    // cx - B^ * u_(n+1/2)

    Teuchos::RCP<Epetra_Vector> cx = DomainExtractor().extract_vector(x, 3);

    Epetra_Vector inter(cx->Map());
    ConStructOp.Multiply(false, *sy, inter);
    cx->Update(-1.0, inter, 1.0);
    ConFluidOp.Multiply(false, *fy, inter);
    cx->Update(-1.0, inter, 1.0);

    Teuchos::RCP<Epetra_Vector> interconsol =
        Teuchos::rcp(new Epetra_Vector(ConStructOp.RangeMap()));
    Core::LinAlg::SolverParams solver_params;
    solver_params.refactor = true;
    solver_params.reset = true;
    constraintsolver_->Solve(interconA_->EpetraOperator(), interconsol, cx, solver_params);

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
    invDiag_->Matrix(0, 0).Multiply(false, *temp1, *temp2);
    sy->Update(-1.0, *temp2, 1.0);

    temp1 = Teuchos::rcp(new Epetra_Vector(fy->Map()));
    temp2 = Teuchos::rcp(new Epetra_Vector(fy->Map()));
    FluidConOp.Multiply(false, *interconsol, *temp1);
    temp1->Scale(alpha_);
    invDiag_->Matrix(1, 1).Multiply(false, *temp1, *temp2);
    fy->Update(-1.0, *temp2, 1.0);
  }

  RangeExtractor().insert_vector(*sy, 0, y);
  RangeExtractor().insert_vector(*fy, 1, y);
  RangeExtractor().insert_vector(*ay, 2, y);
  RangeExtractor().insert_vector(*cy, 3, y);
}



Teuchos::RCP<Core::LinAlg::SparseMatrix> FSI::LungSchurComplement::CalculateSchur(
    const Core::LinAlg::SparseMatrix& A, const Core::LinAlg::SparseMatrix& B,
    const Core::LinAlg::SparseMatrix& C)
{
  // make sure fill_complete was called on the matrices
  if (!A.Filled()) FOUR_C_THROW("A has to be fill_complete");
  if (!B.Filled()) FOUR_C_THROW("B has to be fill_complete");
  if (!C.Filled()) FOUR_C_THROW("C has to be fill_complete");

  temp_ = Core::LinAlg::MLMultiply(A, B, true);
  res_ = Core::LinAlg::MLMultiply(*temp_, C, true);

  return res_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::LungOverlappingBlockMatrix::Label() const
{
  return "FSI::LungOverlappingBlockMatrix";
}

FOUR_C_NAMESPACE_CLOSE
