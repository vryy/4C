/*----------------------------------------------------------------------*/
/*!
 \file poro_monolithicsplit_nopenetration.H

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#ifndef PORO_MONOLITHICSPLIT_NOPENETRATION_H_
#define PORO_MONOLITHICSPLIT_NOPENETRATION_H_


/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/

#include "poro_monolithicsplit.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/

namespace POROELAST
{
  /// monolithic structure split for condensing DOFs, when using the brinkman-equation
  class MonolithicSplitNoPenetration: public MonolithicSplit
  {
  public:
    //! create using a Epetra_Comm
    explicit MonolithicSplitNoPenetration(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! virtual destructor to support polymorph destruction
    virtual ~MonolithicSplitNoPenetration(){};

    //! Setup the monolithic system
    virtual void SetupSystem();

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(bool firstcall=false);

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat);

  protected:
    //! @name Apply current field state to system

    //! Evaluate mechanical-fluid system matrix
    virtual void ApplyStrCouplMatrix(Teuchos::RCP<LINALG::SparseOperator> k_sf //!< mechanical-fluid stiffness matrix
        );

    //! Evaluate fluid-mechanical system matrix
    virtual void ApplyFluidCouplMatrix(Teuchos::RCP<LINALG::SparseOperator> k_fs //!< fluid-mechanical tangent matrix
          );

    //@}
  private:

    //! build block vector from field vectors, e.g. rhs, increment vector
    void SetupVector(Epetra_Vector &f, //!< vector of length of all dofs
        Teuchos::RCP<const Epetra_Vector> sv, //!< vector containing only structural dofs
        Teuchos::RCP<const Epetra_Vector> fv //!< vector containing only fluid dofs
        );

    /// recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time step
    /// (i.e. condensed forces onto the structure) needed for rhs in next time step
    virtual void RecoverLagrangeMultiplier();

    //! @name Global matrixes and vectors

    Teuchos::RCP<LINALG::SparseOperator> k_struct_;
    Teuchos::RCP<LINALG::SparseOperator> k_fluid_;

    Teuchos::RCP<LINALG::SparseOperator> k_lambda_;
    Teuchos::RCP<LINALG::SparseOperator> k_D_;

    Teuchos::RCP<Epetra_Vector> nopenetration_rhs_;

    //@}
  };

}


#endif /* PORO_MONOLITHICSPLIT_NOPENETRATION_H_ */
