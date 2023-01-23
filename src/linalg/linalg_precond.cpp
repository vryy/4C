/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation

\level 0


*----------------------------------------------------------------------*/


#include <Ifpack.h>
#include <Teuchos_TimeMonitor.hpp>

#include "linalg_precond.H"
#include "solver_linalg_solver.H"
#include "linalg_utils_sparse_algebra_manipulation.H"
#include "linalg_mlapi_operator.H"

#include "lib_node.H"
#include "lib_discret.H"



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::Preconditioner::Preconditioner(Teuchos::RCP<Solver> solver) : solver_(solver), ncall_(0) {}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Preconditioner::Setup(Teuchos::RCP<Epetra_Operator> matrix,
    Teuchos::RCP<LINALG::MapExtractor> fsidofmapex, Teuchos::RCP<DRT::Discretization> fdis,
    Teuchos::RCP<Epetra_Map> inodes, bool structuresplit)
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::Preconditioner::Setup");

  Teuchos::Time timer("", true);
  timer.reset();

  std::string solvertype = solver_->Params().get("solver", "none");
  if (solvertype == "aztec" || solvertype == "belos")
  {
    Teuchos::ParameterList* azlist_ptr = NULL;
    if (solvertype == "aztec")
      azlist_ptr = &(solver_->Params().sublist("Aztec Parameters"));
    else
      azlist_ptr = &(solver_->Params().sublist("Belos Parameters"));
    Teuchos::ParameterList& azlist = *azlist_ptr;
    // Teuchos::ParameterList& azlist = solver_->Params().sublist("Aztec Parameters");
    // see whether Operator is a Epetra_CrsMatrix
    Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(&*matrix);

    // get type of preconditioner and build either Ifpack or ML
    // if we have an ifpack parameter list, we do ifpack
    // if we have an ml parameter list we do ml
    bool doifpack = solver_->Params().isSublist("IFPACK Parameters");
    bool doml = solver_->Params().isSublist("ML Parameters");

    if (!A)
    {
      doifpack = false;
      doml = false;
    }

    if (doifpack == false && doml == false)
    {
      dserror("You have to use either ML or Ifpack. No ML Parameters of IFPACK Parameters found!");
    }

    // do ifpack if desired
    if (doifpack)
    {
      prec_ = Teuchos::null;
      Teuchos::ParameterList& ifpacklist = solver_->Params().sublist("IFPACK Parameters");
      ifpacklist.set<bool>("relaxation: zero starting solution", true);
      // create a copy of the scaled matrix
      // so we can reuse the preconditioner
      Pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));
      // get the type of ifpack preconditioner from aztec
      std::string prectype = azlist.get("preconditioner", "ILU");
      int overlap = azlist.get("AZ_overlap", 0);
      Ifpack Factory;
      Ifpack_Preconditioner* prec = Factory.Create(prectype, Pmatrix_.get(), overlap);
      prec->SetParameters(ifpacklist);
      prec->Initialize();
      prec->Compute();
      prec_ = Teuchos::rcp(prec);
    }

    // do ml if desired
    if (doml)
    {
      Teuchos::ParameterList& mllist = solver_->Params().sublist("ML Parameters");
      // see whether we use standard ml or our own mlapi operator
      const bool domlapioperator = mllist.get<bool>("LINALG::AMG_Operator", false);
      if (domlapioperator)
      {
        // create a copy of the scaled matrix
        // so we can reuse the preconditioner several times
        prec_ = Teuchos::null;
        Pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));
        prec_ = Teuchos::rcp(new LINALG::AMG_Operator(Pmatrix_, mllist, true));
      }
      else
      {
        // create a copy of the scaled (and downwinded) matrix
        // so we can reuse the preconditioner several times
        prec_ = Teuchos::null;
        Pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));
        prec_ = Teuchos::rcp(new ML_Epetra::MultiLevelPreconditioner(*Pmatrix_, mllist, true));
        // for debugging ML
        // dynamic_cast<ML_Epetra::MultiLevelPreconditioner&>(*P_).PrintUnused(0);
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Preconditioner::Solve(Teuchos::RCP<Epetra_Operator> matrix,
    Teuchos::RCP<Epetra_MultiVector> x, Teuchos::RCP<Epetra_MultiVector> b, bool refactor,
    bool reset)
{
  std::string solvertype = solver_->Params().get("solver", "none");
  if (solvertype == "aztec" or solvertype == "belos")
  {
    // do just the preconditioner from iterative solver

    // apply the preconditioner
    // This is were the work happens.
    ApplyInverse(*b, *x);
  }
  else
  {
    // this is to bypass an amesos bug that demands the rhs and solution to be
    // the SAME vector in every reuse of the factorization
    // as we can not guarantee that x and b are always the physically same vector,
    // they are always copied to x_ and b_ when the factorization is reused.
    if (refactor || reset)
    {
      b_ = Teuchos::rcp(new Epetra_MultiVector(*b));
      x_ = Teuchos::rcp(new Epetra_MultiVector(*x));
    }
    else
    {
      b_->Update(1.0, *b, 0.0);
      x_->Update(1.0, *x, 0.0);
    }
    // direct solves are done by the solver itself.
    solver_->Solve(matrix, x_, b_, refactor, reset);
    x->Update(1.0, *x_, 0.0);
  }

  ncall_ += 1;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::ParameterList& LINALG::Preconditioner::Params() { return solver_->Params(); }



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::Preconditioner::SetUseTranspose(bool UseTranspose)
{
  return prec_->SetUseTranspose(UseTranspose);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::Preconditioner::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return prec_->Apply(X, Y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::Preconditioner::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return prec_->ApplyInverse(X, Y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double LINALG::Preconditioner::NormInf() const { return prec_->NormInf(); }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* LINALG::Preconditioner::Label() const { return "LINALG::Preconditioner"; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool LINALG::Preconditioner::UseTranspose() const { return prec_->UseTranspose(); }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool LINALG::Preconditioner::HasNormInf() const { return prec_->HasNormInf(); }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Comm& LINALG::Preconditioner::Comm() const { return prec_->Comm(); }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& LINALG::Preconditioner::OperatorDomainMap() const
{
  return prec_->OperatorDomainMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& LINALG::Preconditioner::OperatorRangeMap() const
{
  return prec_->OperatorRangeMap();
}
