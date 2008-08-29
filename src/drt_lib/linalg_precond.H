/*!----------------------------------------------------------------------
\file linalg_precond.H

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef LINALG_PRECOND_H
#define LINALG_PRECOND_H

#include "linalg_sparsematrix.H"

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <stdio.h>


namespace LINALG
{

  class Solver;


  /// linear preconditioning operator
  /*!

    We need preconditioning on matrix blocks for block systems.

    Block systems (for example the monolithic FSI matrix) have to be solved
    using Krylov methods, e.g. GMRES. Those Krylov methods need (block)
    preconditioners. However it is not possible to use Krylov solver on those
    matrix blocks inside the preconditioner, because the outer Krylov method
    demands a linear operator for its preconditioner. Thus we need a way to
    apply just the preconditioners (or a stationary Richardson iteration, a
    direct solver or the like) to our block matrix. This is what this class
    provides.

    \note This class is not used by LINALG::Solver. We use the same input
    parameters here, but that is all. No connection.

    \author u.kue
    \date 05/08
   */
  class Preconditioner : public Epetra_Operator
  {
  public:

    /// construction from a solver object
    explicit Preconditioner(Teuchos::RCP<Solver> solver);

    /// create internal preconditioner object
    /*!
      Destroy any existing one.
     */
    void Setup(Teuchos::RCP<Epetra_Operator> matrix);

    /// Solve system of equations
    /*!
      \param matrix (in/out): system of equations
      \param x      (in/out): initial guess on input, solution on output
      \param b      (in)    : right hand side vector
      \param refactor (in)  : flag indicating whether system should be refactorized
      \param reset  (in)    : flag indicating whether all data from previous solves should
                              be recalculated including preconditioners
     */
    void Solve(Teuchos::RCP<Epetra_Operator>  matrix,
               Teuchos::RCP<Epetra_Vector>    x,
               Teuchos::RCP<Epetra_Vector>    b,
               bool refactor,
               bool reset = false);

    /// @name Attribute set methods

    /// If set true, transpose of this operator will be applied.
    virtual int SetUseTranspose(bool UseTranspose);

    /// @name Mathematical functions

    /// Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    virtual int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    virtual int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

    /// Returns the infinity norm of the global matrix.
    virtual double NormInf() const;

    /// @name Attribute access functions

    /// Returns a character string describing the operator.
    virtual const char* Label() const;

    /// Returns the current UseTranspose setting.
    virtual bool UseTranspose() const;

    /// Returns true if the this object can provide an approximate Inf-norm, false otherwise.
    virtual bool HasNormInf() const;

    /// Returns a pointer to the Epetra_Comm communicator associated with this operator.
    virtual const Epetra_Comm& Comm() const;

    /// Returns the Epetra_Map object associated with the domain of this operator.
    virtual const Epetra_Map& OperatorDomainMap() const;

    /// Returns the Epetra_Map object associated with the range of this operator.
    virtual const Epetra_Map& OperatorRangeMap() const;

    ///@}

  private:

    /// Get number of solver calls done on this solver
    inline int Ncall() const { return ncall_; }

    /// my internal preconditioner
    Teuchos::RCP<Epetra_Operator> prec_;

    //! system of equations used for preconditioning used by prec_ only
    Teuchos::RCP<Epetra_RowMatrix> Pmatrix_;

    /// there is always a solver object
    Teuchos::RCP<Solver> solver_;

    //! counting how many times matrix was solved between resets
    int ncall_;
  };
}

#endif
#endif
