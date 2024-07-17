/*----------------------------------------------------------------------*/
/*! \file

\level 1


\brief Preconditioner for FSI problems with additional constraints
*/
/*----------------------------------------------------------------------*/

#include "4C_fsi_constr_overlapprec.hpp"

#include "4C_adapter_fld_fluid.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_fsi_debugwriter.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_preconditioner_linalg.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

// /*----------------------------------------------------------------------*
//  *----------------------------------------------------------------------*/
FSI::ConstrOverlappingBlockMatrix::ConstrOverlappingBlockMatrix(
    const Core::LinAlg::MultiMapExtractor& maps, Adapter::FSIStructureWrapper& structure,
    Adapter::Fluid& fluid, Adapter::AleFsiWrapper& ale, bool structuresplit, int symmetric,
    double omega, int iterations, double somega, int siterations, double fomega, int fiterations,
    double aomega, int aiterations)
    : FSI::OverlappingBlockMatrix(Teuchos::null, maps, structure, fluid, ale, structuresplit,
          symmetric, omega, iterations, somega, siterations, fomega, fiterations, aomega,
          aiterations)
{
  // determine map of all dofs not related to constraint

  std::vector<Teuchos::RCP<const Epetra_Map>> fsimaps;
  fsimaps.push_back(maps.Map(0));
  fsimaps.push_back(maps.Map(1));
  fsimaps.push_back(maps.Map(2));
  overallfsimap_ = Core::LinAlg::MultiMapExtractor::merge_maps(fsimaps);
  fsiextractor_ = Core::LinAlg::MultiMapExtractor(*overallfsimap_, fsimaps);

  // stuff needed for SIMPLE preconditioner -> this needs to be read
  // in from the input file one day!
  const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();
  alpha_ = fsidyn.sublist("CONSTRAINT").get<double>("ALPHA");
  simpleiter_ = fsidyn.sublist("CONSTRAINT").get<int>("SIMPLEITER");
  prec_ = Core::UTILS::IntegralValue<Inpar::FSI::PrecConstr>(
      fsidyn.sublist("CONSTRAINT"), "PRECONDITIONER");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::ConstrOverlappingBlockMatrix::setup_preconditioner()
{
#ifdef BLOCKMATRIXMERGE

  // this is really evil
  // note that this only works with UMFPACK as fluid solver (saddle
  // point problem!)

  sparse_ = Merge();
  fluidsolver_->setup(sparse_->EpetraMatrix());

#else

  const Core::LinAlg::SparseMatrix& structInnerOp = matrix(0, 0);
  const Core::LinAlg::SparseMatrix& fluidInnerOp = matrix(1, 1);
  const Core::LinAlg::SparseMatrix& aleInnerOp = matrix(2, 2);

  Teuchos::RCP<Core::LinAlg::MapExtractor> fsidofmapex = Teuchos::null;
  Teuchos::RCP<Epetra_Map> irownodes = Teuchos::null;

  structuresolver_->setup(structInnerOp.epetra_matrix());
  fluidsolver_->setup(fluidInnerOp.epetra_matrix(), fsidofmapex, fluid_.discretization(), irownodes,
      structuresplit_);
  if (constalesolver_ == Teuchos::null) alesolver_->setup(aleInnerOp.epetra_matrix());

#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::ConstrOverlappingBlockMatrix::sgs(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  const Core::LinAlg::SparseMatrix& StructInnerOp = matrix(0, 0);
  const Core::LinAlg::SparseMatrix& StructBoundOp = matrix(0, 1);
  const Core::LinAlg::SparseMatrix& StructConOp = matrix(0, 3);
  const Core::LinAlg::SparseMatrix& FluidBoundOp = matrix(1, 0);
  const Core::LinAlg::SparseMatrix& FluidInnerOp = matrix(1, 1);
  const Core::LinAlg::SparseMatrix& FluidMeshOp = matrix(1, 2);
  const Core::LinAlg::SparseMatrix& FluidConOp = matrix(1, 3);
  const Core::LinAlg::SparseMatrix& AleInnerOp = matrix(2, 2);
  const Core::LinAlg::SparseMatrix& ConStructOp = matrix(3, 0);
  const Core::LinAlg::SparseMatrix& ConFluidOp = matrix(3, 1);


  // Extract vector blocks
  // RHS

  const Epetra_Vector& x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // initial guess

  Epetra_Vector& y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> sy = range_extractor().extract_vector(y, 0);
  Teuchos::RCP<Epetra_Vector> fy = range_extractor().extract_vector(y, 1);
  Teuchos::RCP<Epetra_Vector> ay = range_extractor().extract_vector(y, 2);
  Teuchos::RCP<Epetra_Vector> cy = range_extractor().extract_vector(y, 3);

  Teuchos::RCP<Epetra_Vector> sz = Teuchos::rcp(new Epetra_Vector(sy->Map()));
  Teuchos::RCP<Epetra_Vector> fz = Teuchos::rcp(new Epetra_Vector(fy->Map()));
  Teuchos::RCP<Epetra_Vector> az = Teuchos::rcp(new Epetra_Vector(ay->Map()));

  Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp(new Epetra_Vector(domain_map(0)));
  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(domain_map(1)));
  Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(domain_map(2)));

  for (int outerrun = 0; outerrun < simpleiter_; ++outerrun)
  {
    // -------------------------------------------------------------------
    // intermediate fsi dofs: u_(n+1/2) = F^(-1)(f-B^T*lambda_n)
    // -------------------------------------------------------------------

    // inner Richardson loop (FSI block)

    for (int run = 0; run < iterations_; ++run)
    {
      Teuchos::RCP<Epetra_Vector> sx = domain_extractor().extract_vector(x, 0);
      Teuchos::RCP<Epetra_Vector> fx = domain_extractor().extract_vector(x, 1);
      Teuchos::RCP<Epetra_Vector> ax = domain_extractor().extract_vector(x, 2);

      // ----------------------------------------------------------------
      // lower GS


      // Structure
      {
        if (run > 0 or outerrun > 0)
        {
          StructConOp.multiply(false, *cy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);

          StructInnerOp.multiply(false, *sy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);
          StructBoundOp.multiply(false, *fy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);
        }

        // Solve structure equations for sy with the rhs sx
        structuresolver_->solve(StructInnerOp.epetra_matrix(), sz, sx, true);

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
          AleInnerOp.multiply(false, *ay, *tmpax);
          ax->Update(-1.0, *tmpax, 1.0);
        }
        if (structuresplit_)
        {
          if (run > 0 or outerrun > 0)
          {
            const Core::LinAlg::SparseMatrix& aleBoundOp = matrix(2, 1);
            aleBoundOp.multiply(false, *fy, *tmpax);
            ax->Update(-1.0, *tmpax, 1.0);
          }
        }
        else
        {
          const Core::LinAlg::SparseMatrix& aleBoundOp = matrix(2, 0);
          aleBoundOp.multiply(false, *sy, *tmpax);
          ax->Update(-1.0, *tmpax, 1.0);
        }

        alesolver_->solve(AleInnerOp.epetra_matrix(), az, ax, true);

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
          FluidInnerOp.multiply(false, *fy, *tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);

          FluidConOp.multiply(false, *cy, *tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);
        }

        FluidBoundOp.multiply(false, *sy, *tmpfx);
        fx->Update(-1.0, *tmpfx, 1.0);
        FluidMeshOp.multiply(false, *ay, *tmpfx);
        fx->Update(-1.0, *tmpfx, 1.0);


        fluidsolver_->solve(FluidInnerOp.epetra_matrix(), fz, fx, true);

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

      if (symmetric_)
      {
        sx = domain_extractor().extract_vector(x, 0);
        fx = domain_extractor().extract_vector(x, 1);
        ax = domain_extractor().extract_vector(x, 2);

        // ----------------------------------------------------------------
        // upper GS

        {
          FluidInnerOp.multiply(false, *fy, *tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);
          FluidBoundOp.multiply(false, *sy, *tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);
          FluidMeshOp.multiply(false, *ay, *tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);
          FluidConOp.multiply(false, *cy, *tmpfx);
          fx->Update(-1.0, *tmpfx, 1.0);

          fluidsolver_->solve(FluidInnerOp.epetra_matrix(), fz, fx, true);

          local_block_richardson(
              fluidsolver_, FluidInnerOp, fx, fz, tmpfx, fiterations_, fomega_, err_, Comm());
          fy->Update(omega_, *fz, 1.0);
        }

        {
          AleInnerOp.multiply(false, *ay, *tmpax);
          ax->Update(-1.0, *tmpax, 1.0);

          if (structuresplit_)
          {
            const Core::LinAlg::SparseMatrix& aleBoundOp = matrix(2, 1);
            aleBoundOp.multiply(false, *fy, *tmpax);
            ax->Update(-1.0, *tmpax, 1.0);
          }
          else
          {
            const Core::LinAlg::SparseMatrix& aleBoundOp = matrix(2, 0);
            aleBoundOp.multiply(false, *sy, *tmpax);
            ax->Update(-1.0, *tmpax, 1.0);
          }

          alesolver_->solve(AleInnerOp.epetra_matrix(), az, ax, true);

          local_block_richardson(
              alesolver_, AleInnerOp, ax, az, tmpax, aiterations_, aomega_, err_, Comm());

          ay->Update(omega_, *az, 1.0);
        }

        {
          StructInnerOp.multiply(false, *sy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);
          StructBoundOp.multiply(false, *fy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);
          StructConOp.multiply(false, *cy, *tmpsx);
          sx->Update(-1.0, *tmpsx, 1.0);


          structuresolver_->solve(StructInnerOp.epetra_matrix(), sz, sx, true);

          local_block_richardson(
              structuresolver_, StructInnerOp, sx, sz, tmpsx, siterations_, somega_, err_, Comm());
          sy->Update(omega_, *sz, 1.0);
        }
      }
    }

    // -----------------------------------------------------------------------
    // intermediate constraint dofs: Dlambda~ = S^(-1) * (cx - B^ * u_(n+1/2))
    // -----------------------------------------------------------------------

    // cx - B^ * u_(n+1/2)

    Teuchos::RCP<Epetra_Vector> cx = domain_extractor().extract_vector(x, 3);

    Epetra_Vector inter(cx->Map());
    ConStructOp.multiply(false, *sy, inter);
    cx->Update(-1.0, inter, 1.0);
    ConFluidOp.multiply(false, *fy, inter);
    cx->Update(-1.0, inter, 1.0);


    Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy> invDiag(
        fsiextractor_, fsiextractor_, 1, false, true);

    if (prec_ == Inpar::FSI::Simple)
    {
      // D^{-1} = diag(A(0,0))^{-1}

      Epetra_Vector structDiagVec(StructInnerOp.row_map(), false);
      StructInnerOp.extract_diagonal_copy(structDiagVec);
      int err = structDiagVec.Reciprocal(structDiagVec);
      if (err) FOUR_C_THROW("Epetra_MultiVector::Reciprocal (structure matrix) returned %d", err);
      Core::LinAlg::SparseMatrix invstructDiag(structDiagVec);

      Epetra_Vector fluidDiagVec(FluidInnerOp.row_map(), false);
      FluidInnerOp.extract_diagonal_copy(fluidDiagVec);
      err = fluidDiagVec.Reciprocal(fluidDiagVec);
      if (err) FOUR_C_THROW("Epetra_MultiVector::Reciprocal (fluid matrix) returned %d", err);
      Core::LinAlg::SparseMatrix invfluidDiag(fluidDiagVec);

      Epetra_Vector aleDiagVec(AleInnerOp.row_map(), false);
      AleInnerOp.extract_diagonal_copy(aleDiagVec);
      err = aleDiagVec.Reciprocal(aleDiagVec);
      if (err) FOUR_C_THROW("Epetra_MultiVector::Reciprocal (ale matrix) returned %d", err);
      Core::LinAlg::SparseMatrix invaleDiag(aleDiagVec);

      invDiag.assign(0, 0, Core::LinAlg::View, invstructDiag);
      invDiag.assign(1, 1, Core::LinAlg::View, invfluidDiag);
      invDiag.assign(2, 2, Core::LinAlg::View, invaleDiag);
    }
    else if (prec_ == Inpar::FSI::Simplec)
    {
      // D^{-1} = sum(abs(A(0,0)))^{-1}

      Epetra_Vector structDiagVec(StructInnerOp.row_map(), false);
      int err = StructInnerOp.epetra_matrix()->InvRowSums(structDiagVec);
      if (err) FOUR_C_THROW("Epetra_CrsMatrix::InvRowSums (structure matrix) returned %d", err);
      Core::LinAlg::SparseMatrix invstructDiag(structDiagVec);

      Epetra_Vector fluidDiagVec(FluidInnerOp.row_map(), false);
      err = FluidInnerOp.epetra_matrix()->InvRowSums(fluidDiagVec);
      if (err) FOUR_C_THROW("Epetra_CrsMatrix::InvRowSums (fluid matrix) returned %d", err);
      Core::LinAlg::SparseMatrix invfluidDiag(fluidDiagVec);

      Epetra_Vector aleDiagVec(AleInnerOp.row_map(), false);
      err = AleInnerOp.epetra_matrix()->InvRowSums(aleDiagVec);
      if (err) FOUR_C_THROW("Epetra_CrsMatrix::InvRowSums (ale matrix) returned %d", err);
      Core::LinAlg::SparseMatrix invaleDiag(aleDiagVec);

      invDiag.assign(0, 0, Core::LinAlg::View, invstructDiag);
      invDiag.assign(1, 1, Core::LinAlg::View, invfluidDiag);
      invDiag.assign(2, 2, Core::LinAlg::View, invaleDiag);
    }
    else
      FOUR_C_THROW("Unknown type of preconditioner for constraint fsi system");

    invDiag.complete();


    // S = - B^ * D^{-1} * B^T

    Teuchos::RCP<Core::LinAlg::SparseMatrix> temps =
        Core::LinAlg::MatrixMultiply(ConStructOp, false, invDiag.matrix(0, 0), false);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> interconA =
        Core::LinAlg::MatrixMultiply(*temps, false, StructConOp, false, false);

    temps = Core::LinAlg::MatrixMultiply(ConFluidOp, false, invDiag.matrix(1, 1), false);
    Teuchos::RCP<Core::LinAlg::SparseMatrix> tempss =
        Core::LinAlg::MatrixMultiply(*temps, false, FluidConOp, false);
    interconA->add(*tempss, false, 1.0, 1.0);
    interconA->complete(StructConOp.domain_map(), ConStructOp.range_map());
    interconA->scale(-1.0);

    Teuchos::ParameterList constrsolvparams;
    Core::UTILS::AddEnumClassToParameterList<Core::LinearSolver::SolverType>(
        "SOLVER", Core::LinearSolver::SolverType::umfpack, constrsolvparams);
    Teuchos::RCP<Epetra_Vector> interconsol =
        Teuchos::rcp(new Epetra_Vector(ConStructOp.range_map()));
    Teuchos::RCP<Core::LinAlg::Solver> ConstraintSolver = Teuchos::rcp(new Core::LinAlg::Solver(
        constrsolvparams, interconA->Comm(), Global::Problem::instance()->solver_params_callback(),
        Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
            Global::Problem::instance()->io_params(), "VERBOSITY")));
    Core::LinAlg::SolverParams solver_params;
    solver_params.refactor = true;
    solver_params.reset = true;
    ConstraintSolver->solve(interconA->epetra_operator(), interconsol, cx, solver_params);

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
    StructConOp.multiply(false, *interconsol, *temp1);
    temp1->Scale(alpha_);
    invDiag.matrix(0, 0).multiply(false, *temp1, *temp2);
    sy->Update(-1.0, *temp2, 1.0);

    temp1 = Teuchos::rcp(new Epetra_Vector(fy->Map()));
    temp2 = Teuchos::rcp(new Epetra_Vector(fy->Map()));
    FluidConOp.multiply(false, *interconsol, *temp1);
    temp1->Scale(alpha_);
    invDiag.matrix(1, 1).multiply(false, *temp1, *temp2);
    fy->Update(-1.0, *temp2, 1.0);
  }

  range_extractor().insert_vector(*sy, 0, y);
  range_extractor().insert_vector(*fy, 1, y);
  range_extractor().insert_vector(*ay, 2, y);
  range_extractor().insert_vector(*cy, 3, y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::ConstrOverlappingBlockMatrix::Label() const
{
  return "FSI::ConstrOverlappingBlockMatrix";
}

FOUR_C_NAMESPACE_CLOSE
