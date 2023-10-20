/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 0

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_LINEAR_SOLVER_PRECONDITIONER_LINALG_H
#define BACI_LINEAR_SOLVER_PRECONDITIONER_LINALG_H

#include "baci_linalg_sparsematrix.H"

#include <stdio.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

namespace DRT
{
  class Discretization;
}

namespace CORE::LINALG
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

    \note This class is not used by CORE::LINALG::Solver. We use the same input
    parameters here, but that is all. No connection.

    \author u.kue
    \date 05/08
   */
  class Preconditioner : public Epetra_Operator
  {
   public:
    /// construction from a solver object
    explicit Preconditioner(Teuchos::RCP<Solver> solver);

    ~Preconditioner() override
    {
      // destruction order is important
      solver_ = Teuchos::null;
      prec_ = Teuchos::null;
      Pmatrix_ = Teuchos::null;
    }

    /// create internal preconditioner object
    /*!
      Destroy any existing one.
     */
    void Setup(Teuchos::RCP<Epetra_Operator> matrix,
        Teuchos::RCP<CORE::LINALG::MapExtractor> fsidofmapex = Teuchos::null,
        Teuchos::RCP<::DRT::Discretization> fdis = Teuchos::null,
        Teuchos::RCP<Epetra_Map> inodes = Teuchos::null, bool structuresplit = false);

    /// Solve system of equations
    /*!
      \param matrix (in/out): system of equations
      \param x      (in/out): initial guess on input, solution on output
      \param b      (in)    : right hand side vector
      \param refactor (in)  : flag indicating whether system should be refactorized
      \param reset  (in)    : flag indicating whether all data from previous solves should
                              be recalculated including preconditioners
     */
    void Solve(Teuchos::RCP<Epetra_Operator> matrix, Teuchos::RCP<Epetra_MultiVector> x,
        Teuchos::RCP<Epetra_MultiVector> b, bool refactor, bool reset = false);


    /// get flag from solver whether factorization has been performed
    //     bool IsFactored() const;

    /// get underlying preconditioner Epetra_Operator
    Teuchos::RCP<Epetra_Operator> EpetraOperator() { return prec_; }

    /// get underlying solver parameter list
    Teuchos::ParameterList& Params();

    /// @name Attribute set methods

    /// If set true, transpose of this operator will be applied.
    int SetUseTranspose(bool UseTranspose) override;

    /// @name Mathematical functions

    /// Returns the result of a Epetra_Operator applied to a Epetra_MultiVector X in Y.
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

    /// Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const override;

    /// Returns the infinity norm of the global matrix.
    double NormInf() const override;

    /// @name Attribute access functions

    /// Returns a character string describing the operator.
    const char* Label() const override;

    /// Returns the current UseTranspose setting.
    bool UseTranspose() const override;

    /// Returns true if the this object can provide an approximate Inf-norm, false otherwise.
    bool HasNormInf() const override;

    /// Returns a pointer to the Epetra_Comm communicator associated with this operator.
    const Epetra_Comm& Comm() const override;

    /// Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map& OperatorDomainMap() const override;

    /// Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map& OperatorRangeMap() const override;

    ///@}

   private:
    /// Get number of solver calls done on this solver
    inline int Ncall() const { return ncall_; }

    /// my internal preconditioner
    Teuchos::RCP<Epetra_Operator> prec_;

    //! system of equations used for preconditioning used by prec_ only
    Teuchos::RCP<Epetra_RowMatrix> Pmatrix_;

    /// there is always a solver object
    Teuchos::RCP<Epetra_MultiVector> x_;
    Teuchos::RCP<Epetra_MultiVector> b_;
    Teuchos::RCP<Solver> solver_;

    //! counting how many times matrix was solved between resets
    int ncall_;
  };
}  // namespace CORE::LINALG

#endif
