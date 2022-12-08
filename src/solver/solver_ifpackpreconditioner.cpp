/*----------------------------------------------------------------------------*/
/*! \file

\brief Declaration
\level 0

\brief LINALG::SOLVER wrapper around Trilinos' IFPACK preconditioner
*/
/*----------------------------------------------------------------------------*/

#include "dserror.H"

#include "solver_ifpackpreconditioner.H"

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
LINALG::SOLVER::IFPACKPreconditioner::IFPACKPreconditioner(
    FILE* outfile, Teuchos::ParameterList& ifpacklist, Teuchos::ParameterList& solverlist)
    : PreconditionerType(outfile), ifpacklist_(ifpacklist), solverlist_(solverlist)
{
  return;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void LINALG::SOLVER::IFPACKPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  SetupLinearProblem(matrix, x, b);

  if (create)
  {
    Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(matrix);
    if (A == NULL) dserror("CrsMatrix expected");

    // free old matrix first
    prec_ = Teuchos::null;
    Pmatrix_ = Teuchos::null;

    // create a copy of the scaled matrix
    // so we can reuse the preconditioner
    Pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));

    // get the type of ifpack preconditioner from solver parameter list
    std::string prectype = solverlist_.get("Preconditioner Type", "ILU");
    const int overlap = ifpacklist_.get("IFPACKOVERLAP", 0);

    // create the preconditioner
    Ifpack Factory;
    prec_ = Teuchos::rcp(Factory.Create(prectype, Pmatrix_.get(), overlap));

    if (prec_.is_null())
      dserror("Creation of IFPACK preconditioner of type '%s' failed.", prectype.c_str());

    // setup
    prec_->SetParameters(ifpacklist_);
    prec_->Initialize();
    prec_->Compute();

    return;
  }
}
