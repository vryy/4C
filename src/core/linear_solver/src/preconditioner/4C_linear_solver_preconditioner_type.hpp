/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_LINEAR_SOLVER_PRECONDITIONER_TYPE_HPP
#define FOUR_C_LINEAR_SOLVER_PRECONDITIONER_TYPE_HPP

#include "4C_config.hpp"

#include <Epetra_Comm.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINEAR_SOLVER
{
  /// preconditioner base class
  /*!
     A KrylovSolver object needs one (or more) preconditioner objects. There
     are many possible preconditioners. A unified framework simplifies the
     solution process.
  */
  class PreconditionerTypeBase
  {
   public:
    /*!
       No setup is done upon construction, only the preconditioner object is
       created.
    */
    PreconditionerTypeBase() = default;

    /// virtual destruction
    virtual ~PreconditionerTypeBase() = default;

    /// Setup preconditioner with a given linear system.
    virtual void Setup(
        bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b) = 0;

    /// linear operator used for preconditioning
    virtual Teuchos::RCP<Epetra_Operator> PrecOperator() const = 0;
  };
}  // namespace CORE::LINEAR_SOLVER

FOUR_C_NAMESPACE_CLOSE

#endif
