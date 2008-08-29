/*!----------------------------------------------------------------------
\file linalg_precond.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "linalg_precond.H"
#include "linalg_solver.H"

#include "linalg_mlapi_operator.H"
#include "simpler_operator.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::Preconditioner::Preconditioner(Teuchos::RCP<Solver> solver)
  : solver_(solver),
    ncall_(0)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Preconditioner::Setup(Teuchos::RCP<Epetra_Operator> matrix)
{
  std::string solvertype = solver_->Params().get("solver","none");
  if (solvertype=="aztec")
  {
    Teuchos::ParameterList& azlist = solver_->Params().sublist("Aztec Parameters");
    // see whether Operator is a Epetra_CrsMatrix
    Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(&*matrix);

    // get type of preconditioner and build either Ifpack or ML
    // if we have an ifpack parameter list, we do ifpack
    // if we have an ml parameter list we do ml
    bool   doifpack  = solver_->Params().isSublist("IFPACK Parameters");
    bool   doml      = solver_->Params().isSublist("ML Parameters");
    bool   dosimpler = solver_->Params().isSublist("SIMPLER");
    if (!A || dosimpler)
    {
      doifpack = false;
      doml     = false;
    }

    // do ifpack if desired
    if (doifpack)
    {
      Teuchos::ParameterList& ifpacklist = solver_->Params().sublist("IFPACK Parameters");
      ifpacklist.set<bool>("relaxation: zero starting solution",true);
      // create a copy of the scaled matrix
      // so we can reuse the preconditioner
      Pmatrix_ = rcp(new Epetra_CrsMatrix(*A));
      // get the type of ifpack preconditioner from aztec
      std::string prectype = azlist.get("preconditioner","ILU");
      int    overlap  = azlist.get("AZ_overlap",0);
      Ifpack Factory;
      Ifpack_Preconditioner* prec = Factory.Create(prectype,Pmatrix_.get(),overlap);
      prec->SetParameters(ifpacklist);
      prec->Initialize();
      prec->Compute();
      prec_ = rcp(prec);
    }

    // do ml if desired
    if (doml)
    {
      Teuchos::ParameterList& mllist = solver_->Params().sublist("ML Parameters");
      // see whether we use standard ml or our own mlapi operator
      const bool domlapioperator = mllist.get<bool>("LINALG::AMG_Operator",false);
      if (domlapioperator)
      {
        // create a copy of the scaled matrix
        // so we can reuse the preconditioner several times
        Pmatrix_ = rcp(new Epetra_CrsMatrix(*A));
        prec_ = rcp(new LINALG::AMG_Operator(Pmatrix_,mllist,true));
      }
      else
      {
        // create a copy of the scaled (and downwinded) matrix
        // so we can reuse the preconditioner several times
        Pmatrix_ = rcp(new Epetra_CrsMatrix(*A));
        prec_ = rcp(new ML_Epetra::MultiLevelPreconditioner(*Pmatrix_,mllist,true));
        // for debugging ML
        //dynamic_cast<ML_Epetra::MultiLevelPreconditioner&>(*P_).PrintUnused(0);
      }
    }

#if 0
    if (dosimpler)
    {
      // SIMPLER does not need copy of preconditioning matrix to live
      // SIMPLER does not use the downwinding installed here, it does
      // its own downwinding inside if desired
      prec_ = rcp(new LINALG::SIMPLER_Operator(matrix,Params(),
                                               solver_->Params().sublist("SIMPLER"),
                                               outfile_));
      Pmatrix_ = null;
    }
#endif
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Preconditioner::Solve(Teuchos::RCP<Epetra_Operator>  matrix,
                                   Teuchos::RCP<Epetra_Vector>    x,
                                   Teuchos::RCP<Epetra_Vector>    b,
                                   bool refactor,
                                   bool reset)
{
  std::string solvertype = solver_->Params().get("solver","none");
  if (solvertype=="aztec")
  {
    // do just the preconditioner from iterative solver

#if 0
    bool   doifpack  = solver_->Params().isSublist("IFPACK Parameters");
    if (doifpack)
    {
      Teuchos::ParameterList& ifpacklist = solver_->Params().sublist("IFPACK Parameters");
      ifpacklist.set<bool>("relaxation: zero starting solution",false);
    }
#endif

    // apply the preconditioner
    // This is were the work happens.
    ApplyInverse(*b,*x);
  }
  else
  {
    // direct solves are done by the solver itself.
    solver_->Solve(matrix,x,b,refactor,reset);
  }

  ncall_ += 1;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::Preconditioner::SetUseTranspose(bool UseTranspose)
{
  return prec_->SetUseTranspose(UseTranspose);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::Preconditioner::Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  return prec_->Apply(X,Y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::Preconditioner::ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
{
  return prec_->ApplyInverse(X,Y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double LINALG::Preconditioner::NormInf() const
{
  return prec_->NormInf();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* LINALG::Preconditioner::Label() const
{
  return "LINALG::Preconditioner";
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool LINALG::Preconditioner::UseTranspose() const
{
  return prec_->UseTranspose();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool LINALG::Preconditioner::HasNormInf() const
{
  return prec_->HasNormInf();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Comm& LINALG::Preconditioner::Comm() const
{
  return prec_->Comm();
}


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


#endif
