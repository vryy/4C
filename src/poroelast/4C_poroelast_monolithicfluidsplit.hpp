/*----------------------------------------------------------------------*/
/*! \file

 \brief  monolithic poroelasticity algorithm with split of fluid degrees of freedom at the interface

\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_POROELAST_MONOLITHICFLUIDSPLIT_HPP
#define FOUR_C_POROELAST_MONOLITHICFLUIDSPLIT_HPP


#include "4C_config.hpp"

#include "4C_poroelast_monolithicsplit.hpp"

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class MatrixRowTransform;
  class MatrixColTransform;
  class MatrixRowColTransform;
}  // namespace CORE::LINALG

namespace POROELAST
{
  //! monolithic fluid split for condensing DOFs, when using the brinkman-equation
  class MonolithicFluidSplit : public MonolithicSplit
  {
   public:
    //! create using a Epetra_Comm
    explicit MonolithicFluidSplit(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams,
        Teuchos::RCP<CORE::LINALG::MapExtractor> porosity_splitter);

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
    void SetupSystem() override;

    //! setup composed right hand side from field solvers
    void SetupRHS(bool firstcall = false) override;

    //! setup composed system matrix from field solvers
    void SetupSystemMatrix(CORE::LINALG::BlockSparseMatrixBase& mat) override;

   private:
    //! build block vector from field vectors
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
        bool firstcall = false) override;

    //! recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time
    //! step (i.e. condensed forces onto the structure) needed for rhs in next time step
    void RecoverLagrangeMultiplierAfterTimeStep() override;

    //! @name matrix transformation
    //! transform object for fluid interface matrix \f$F_{\Gamma \Gamma}\f$
    Teuchos::RCP<CORE::LINALG::MatrixRowColTransform> fggtransform_;
    //! transform object for fluid interface matrix \f$F_{\Gamma I}\f$
    Teuchos::RCP<CORE::LINALG::MatrixRowTransform> fgitransform_;
    //! transform object for fluid interface matrix \f$F_{I \Gamma}\f$
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> figtransform_;
    //! transform object for fluid coupling matrix \f$C_{\Gamma \Gamma}^F\f$
    Teuchos::RCP<CORE::LINALG::MatrixRowTransform> cfggtransform_;
    //! transform object for structure coupling matrix \f$C_{\Gamma \Gamma}^S\f$
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> csggtransform_;
    //! transform object for fluid coupling matrix \f$C_{\Gamma I}^F\f$
    Teuchos::RCP<CORE::LINALG::MatrixRowTransform> cfgitransform_;
    //! transform object for structure coupling matrix \f$C_{I \Gamma}^S\f$
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> csigtransform_;
    //!@}

    //! block \f$F_{\Gamma I,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const CORE::LINALG::SparseOperator> fgicur_;

    //! block \f$F_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const CORE::LINALG::SparseOperator> fggcur_;

    //! block \f$C_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const CORE::LINALG::SparseOperator> cgicur_;

    //! block \f$C_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const CORE::LINALG::SparseOperator> cggcur_;

    //!@}
  };

}  // namespace POROELAST

FOUR_C_NAMESPACE_CLOSE

#endif
