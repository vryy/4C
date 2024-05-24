/*----------------------------------------------------------------------------*/
/*! \file

\brief Declaration
\level 0

\brief CORE::LINALG::SOLVER wrapper around Trilinos' IFPACK preconditioner
*/
/*----------------------------------------------------------------------------*/

#include "4C_linear_solver_preconditioner_ifpack.hpp"

#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
CORE::LINEAR_SOLVER::IFPACKPreconditioner::IFPACKPreconditioner(
    Teuchos::ParameterList& ifpacklist, Teuchos::ParameterList& solverlist)
    : ifpacklist_(ifpacklist), solverlist_(solverlist)
{
  return;
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::IFPACKPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  if (create)
  {
    Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>(matrix);
    if (A == nullptr) FOUR_C_THROW("CrsMatrix expected");

    // free old matrix first
    pmatrix_ = Teuchos::null;

    // create a copy of the scaled matrix
    // so we can reuse the preconditioner
    pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));

    // get the type of ifpack preconditioner from solver parameter list
    std::string prectype = solverlist_.get("Preconditioner Type", "ILU");
    const int overlap = ifpacklist_.get("IFPACKOVERLAP", 0);

    // create the preconditioner
    Ifpack Factory;
    Teuchos::RCP<Ifpack_Preconditioner> prec =
        Teuchos::rcp(Factory.Create(prectype, pmatrix_.get(), overlap));

    if (prec.is_null())
      FOUR_C_THROW("Creation of IFPACK preconditioner of type '%s' failed.", prectype.c_str());

    // setup
    prec->SetParameters(ifpacklist_);
    prec->Initialize();
    prec->Compute();

    preconditioner_operator_ = prec;

    return;
  }
}

FOUR_C_NAMESPACE_CLOSE
