/*!----------------------------------------------------------------------
\file solver_preconditionertype.H

<pre>
\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
Created on: Jul 4, 2011
</pre>
*/
/*----------------------------------------------------------------------*/
#ifndef SOLVER_PRECONDITIONERTYPE_H_
#define SOLVER_PRECONDITIONERTYPE_H_

#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace LINALG
{
  namespace SOLVER
  {
    /// preconditioner base class
    /*!

      A KrylovSolver object needs one (or more) preconditioner objects. There
      are many possible preconditioners. A unified framework simplifies the
      solution process.

      \author u.kue
      \date 04/11
     */
    class PreconditionerType
    {
     public:
      /// construction with output on proc 0 (might be NULL)
      /*!
        No setup is done upon construction, only the preconditioner object is
        created.
       */
      explicit PreconditionerType(FILE* outfile);

      /// virtual destruction
      virtual ~PreconditionerType() {}

      /// linear problem created (managed) by this preconditioner
      /*!
        This is how aztec sees the linear problem that needs to be solved.
       */
      virtual Epetra_LinearProblem& LinearProblem() { return lp_; }

      /*! \brief Support routine for setup
       *
       * Pass the components of the linear system on to the Epetra_LinearProblem.
       *
       * @param matrix Matrix of the linear problem
       * @param x Solution vector of the linear problem
       * @param b Right-hand side of the linear problem
       */
      void SetupLinearProblem(
          Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b);

      /// Setup preconditioner with a given linear system.
      virtual void Setup(
          bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b) = 0;

      /// Finish calculation after linear solve.
      /*!
        This is empty in most cases, however some preconditioners might want to
        scale the solution.
       */
      virtual void Finish(Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b) {}

      /// linear operator used for preconditioning
      /*!
        This is how aztec sees the linear operator to be used for preconditioning.
       */
      virtual Epetra_Operator* PrecOperator() const = 0;

      /// linear operator used for preconditioning
      virtual Teuchos::RCP<Epetra_Operator> PrecOperatorRCP() const = 0;

      /// Apply inverse of the preconditioner
      int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
      {
        return PrecOperator()->ApplyInverse(X, Y);
      }

      /// debug output
      virtual void Print(std::ostream& stream) = 0;

      /// return name of sublist in paramterlist which contains parameters for preconditioner
      /// e.g. "IFPACK parameters" for an Ifpack preconditioner
      virtual const std::string getParameterListName() const { return "unknown"; }

     protected:
      //! a linear problem wrapper class used by Trilinos and for scaling of the system
      Epetra_LinearProblem lp_;

      //! file to write output to (proc 0 only, can be NULL on input)
      FILE* outfile_;
    };
  }  // namespace SOLVER
}  // namespace LINALG


#endif /* SOLVER_PRECONDITIONERTYPE_H_ */
