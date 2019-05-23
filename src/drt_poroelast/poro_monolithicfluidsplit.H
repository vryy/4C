/*----------------------------------------------------------------------*/
/*!
 \file poro_monolithicfluidsplit.H

 \brief  monolithic poroelasticity algorithm with split of fluid degrees of freedom at the interface

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *------------------------------------------------------------------------------------------------*/

#ifndef PORO_MONOLITHICFLUIDSPLIT_H_
#define PORO_MONOLITHICFLUIDSPLIT_H_

/*----------------------------------------------------------------------*
 | headers                                                              |
 *----------------------------------------------------------------------*/

#include "poro_monolithicsplit.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/

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

namespace POROELAST
{
  /// monolithic fluid split for condensing DOFs, when using the brinkman-equation
  class MonolithicFluidSplit : public MonolithicSplit
  {
   public:
    //! create using a Epetra_Comm
    explicit MonolithicFluidSplit(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    //! virtual destructor to support polymorph destruction
    virtual ~MonolithicFluidSplit(){};

    /*! do the setup for the monolithic system


     1.) setup coupling
     2.) get maps for all blocks in the system (and for the whole system as well)
     create combined map
     3.) create system matrix


     \note We want to do this setup after reading the restart information, not
     directly in the constructor. This is necessary since during restart (if
     ReadMesh is called), the dofmaps for the blocks might get invalid.
     */
    //! Setup the monolithic system
    virtual void SetupSystem();

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(bool firstcall = false);

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat);

   private:
    /// build block vector from field vectors
    void SetupVector(Epetra_Vector& f, Teuchos::RCP<const Epetra_Vector> sv,
        Teuchos::RCP<const Epetra_Vector> fv, double fluidscale);

    //! extract the field vectors from a given composed vector
    /*!
     \param x  (i) composed vector that contains all field vectors
     \param sx (o) structural vector (e.g. displacements)
     \param fx (o) fluid vector (e.g. velocities and pressure)
     */
    void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
        bool firstcall = false);

    /// recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time
    /// step (i.e. condensed forces onto the structure) needed for rhs in next time step
    virtual void RecoverLagrangeMultiplierAfterTimeStep();

    /// @name matrix transformation
    //! transform object for fluid interface matrix \f$F_{\Gamma \Gamma}\f$
    Teuchos::RCP<FSI::UTILS::MatrixRowColTransform> fggtransform_;
    //! transform object for fluid interface matrix \f$F_{\Gamma I}\f$
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> fgitransform_;
    //! transform object for fluid interface matrix \f$F_{I \Gamma}\f$
    Teuchos::RCP<FSI::UTILS::MatrixColTransform> figtransform_;
    //! transform object for fluid coupling matrix \f$C_{\Gamma \Gamma}^F\f$
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> cfggtransform_;
    //! transform object for structure coupling matrix \f$C_{\Gamma \Gamma}^S\f$
    Teuchos::RCP<FSI::UTILS::MatrixColTransform> csggtransform_;
    //! transform object for fluid coupling matrix \f$C_{\Gamma I}^F\f$
    Teuchos::RCP<FSI::UTILS::MatrixRowTransform> cfgitransform_;
    //! transform object for structure coupling matrix \f$C_{I \Gamma}^S\f$
    Teuchos::RCP<FSI::UTILS::MatrixColTransform> csigtransform_;
    //@}

    //! block \f$F_{\Gamma I,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseOperator> fgicur_;

    //! block \f$F_{\Gamma I,i}\f$ of fluid matrix at previous iteration \f$i\f$
    // Teuchos::RCP<const LINALG::SparseOperator> fgipre_;

    //! block \f$F_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseOperator> fggcur_;

    //! block \f$F_{\Gamma\Gamma,i}\f$ of fluid matrix at previous iteration \f$i\f$
    // Teuchos::RCP<const LINALG::SparseOperator> fggpre_;

    //! block \f$C_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseOperator> cgicur_;

    //! block \f$C_{\Gamma\Gamma,i}\f$ of fluid matrix at previous iteration \f$i\f$
    // Teuchos::RCP<const LINALG::SparseOperator> cgipre_;

    //! block \f$C_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseOperator> cggcur_;

    //! block \f$C_{\Gamma\Gamma,i}\f$ of fluid matrix at previous iteration \f$i\f$
    // Teuchos::RCP<const LINALG::SparseOperator> cggpre_;

    //@}

  };  // MonolithicFluidSplit

}  // namespace POROELAST

#endif /* PORO_MONOLITHICFLUIDSPLIT_H_ */
