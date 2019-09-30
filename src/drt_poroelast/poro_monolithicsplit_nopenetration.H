/*----------------------------------------------------------------------*/
/*! \file

 \brief porous medium algorithm with matrix split for condensation of
      no-penetration constraint

\level 2

\maintainer Johannes Kremheller
            kremheller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *----------------------------------------------------------------------*/


#ifndef PORO_MONOLITHICSPLIT_NOPENETRATION_H_
#define PORO_MONOLITHICSPLIT_NOPENETRATION_H_


/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/

#include "poro_monolithicsplit.H"

namespace FSI
{
  // class OverlappingBlockMatrix;

  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }  // namespace UTILS
}  // namespace FSI

namespace ADAPTER
{
  class CouplingNonLinMortar;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/

namespace POROELAST
{
  /// monolithic structure split for condensing DOFs, when using the brinkman-equation
  class MonolithicSplitNoPenetration : public MonolithicSplit
  {
   public:
    //! create using a Epetra_Comm
    explicit MonolithicSplitNoPenetration(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! virtual destructor to support polymorph destruction
    virtual ~MonolithicSplitNoPenetration(){};

    //! Setup the monolithic system
    virtual void SetupSystem();

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(bool firstcall = false);

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat);

    //! take current results for converged and save for next time step
    virtual void Update();

    /// read restart data
    virtual void ReadRestart(const int step);

    //! contains text to PrintNewtonIter
    virtual void PrintNewtonIterTextStream(std::ostringstream& oss);

    //! contains header to PrintNewtonIter
    virtual void PrintNewtonIterHeaderStream(std::ostringstream& oss);


   protected:
    //! @name Apply current field state to system

    //! Evaluate mechanical-fluid system matrix
    virtual void ApplyStrCouplMatrix(
        Teuchos::RCP<LINALG::SparseOperator> k_sf  //!< mechanical-fluid stiffness matrix
    );

    //! Evaluate fluid-mechanical system matrix
    virtual void ApplyFluidCouplMatrix(
        Teuchos::RCP<LINALG::SparseOperator> k_fs  //!< fluid-mechanical tangent matrix
    );

    /// recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each
    /// iteration step (i.e. condensed forces onto the structure) needed for rhs in next newton step
    virtual void RecoverLagrangeMultiplierAfterNewtonStep(Teuchos::RCP<const Epetra_Vector> x);

    /// recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time
    /// step (i.e. condensed forces onto the structure) needed for rhs in next time step
    virtual void RecoverLagrangeMultiplierAfterTimeStep();

    //! output
    virtual void Output(bool forced_writerestart = false);

    //! setup of coupling object and system matrices
    virtual void SetupCouplingAndMatrices();

    //! start a new time step
    virtual void PrepareTimeStep();

    //! convergence check for Newton solver
    virtual void BuildConvergenceNorms();

    //@}
   private:
    //! build block vector from field vectors, e.g. rhs, increment vector
    void SetupVector(Epetra_Vector& f,         //!< vector of length of all dofs
        Teuchos::RCP<const Epetra_Vector> sv,  //!< vector containing only structural dofs
        Teuchos::RCP<const Epetra_Vector> fv   //!< vector containing only fluid dofs
    );

    //! @name Global matrices and vectors

    Teuchos::RCP<LINALG::SparseOperator> k_struct_;
    Teuchos::RCP<LINALG::SparseOperator> k_fluid_;

    Teuchos::RCP<LINALG::SparseMatrix> k_lambda_;
    Teuchos::RCP<LINALG::SparseMatrix> k_D_;
    Teuchos::RCP<LINALG::SparseMatrix> k_invD_;

    Teuchos::RCP<LINALG::SparseMatrix> k_Dn_;

    Teuchos::RCP<LINALG::SparseMatrix> k_lambdainvD_;

    Teuchos::RCP<LINALG::BlockSparseMatrixBase> k_porodisp_;
    Teuchos::RCP<LINALG::SparseOperator> k_porofluid_;

    Teuchos::RCP<Epetra_Vector> nopenetration_rhs_;

    //! transform object for k_D matrix \f$D\f$
    Teuchos::RCP<FSI::UTILS::MatrixColTransform> k_D_transform_;
    //! transform object for k_D matrix \f$D\f$
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> k_invD_transform_;

    //! transform object for linearization of k_D matrix \f$D\f$
    Teuchos::RCP<FSI::UTILS::MatrixColTransform> k_DLin_transform_;

    //! Lagrange multiplier \f$\lambda_\Gamma^{n+1}\f$ at the interface (ie condensed forces onto
    //! the structure) evaluated at actual iteration step \f$t_{n+1}\f$ but needed for next
    //! iteration step
    Teuchos::RCP<Epetra_Vector> lambdanp_;

    //@}

    /// @name Some quantities to recover the Langrange multiplier at the end of each iteration step

    //! block \f$F_{\Gamma I,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fgicur_;

    //! block \f$S_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fggcur_;

    //! block \f$Cfs_{\Gamma\Gamma,i+1}\f$ of fs-coupling matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> cfsgicur_;

    //! block \f$Cfs_{\Gamma\Gamma,i+1}\f$ of fs-coupling matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> cfsggcur_;

    Teuchos::RCP<Epetra_Vector> rhs_fgcur_;

    //@}

    //! @name Iterative solution technique

    //!< norm of no-penetration constraint
    double normrhs_nopenetration_;
    //@}

    Teuchos::RCP<ADAPTER::CouplingNonLinMortar> mortar_adapter_;
  };

}  // namespace POROELAST


#endif /* PORO_MONOLITHICSPLIT_NOPENETRATION_H_ */
