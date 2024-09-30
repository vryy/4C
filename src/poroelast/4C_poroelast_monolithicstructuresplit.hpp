/*----------------------------------------------------------------------*/
/*! \file

 \brief  monolithic poroelasticity algorithm with split of structure degrees of freedom at the
interface

\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_POROELAST_MONOLITHICSTRUCTURESPLIT_HPP
#define FOUR_C_POROELAST_MONOLITHICSTRUCTURESPLIT_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter_converter.hpp"
#include "4C_poroelast_monolithicsplit.hpp"

FOUR_C_NAMESPACE_OPEN


namespace PoroElast
{
  //! monolithic structure split for condensing DOFs, when using the brinkman-equation
  class MonolithicStructureSplit : public MonolithicSplit
  {
   public:
    //! create using a Epetra_Comm
    explicit MonolithicStructureSplit(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams,
        Teuchos::RCP<Core::LinAlg::MapExtractor> porosity_splitter);

    /*! do the setup for the monolithic system


     1.) setup coupling
     2.) get maps for all blocks in the system (and for the whole system as well)
     create combined map
     3.) create system matrix


     \note We want to do this setup after reading the restart information, not
     directly in the constructor. This is necessary since during restart (if
     read_mesh is called), the dofmaps for the blocks might get invalid.
     */
    //! Setup the monolithic system
    void setup_system() override;

    //! setup composed right hand side from field solvers
    void setup_rhs(bool firstcall = false) override;

    //! setup composed system matrix from field solvers
    void setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat) override;

   private:
    //! build block vector from field vectors
    void setup_vector(Core::LinAlg::Vector& f, Teuchos::RCP<const Core::LinAlg::Vector> sv,
        Teuchos::RCP<const Core::LinAlg::Vector> fv, double fluidscale);

    //! extract the field vectors from a given composed vector
    /*!
     \param x  (i) composed vector that contains all field vectors
     \param sx (o) structural vector (e.g. displacements)
     \param fx (o) fluid vector (e.g. velocities and pressure)
     */
    void extract_field_vectors(Teuchos::RCP<const Core::LinAlg::Vector> x,
        Teuchos::RCP<const Core::LinAlg::Vector>& sx, Teuchos::RCP<const Core::LinAlg::Vector>& fx,
        bool firstcall = false) override;

    //! recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time
    //! step (i.e. condensed forces onto the structure) needed for rhs in next time step
    void recover_lagrange_multiplier_after_time_step() override;

    //! @name matrix transformation
    //! transform object for structure interface matrix \f$S_{\Gamma \Gamma}\f$
    Teuchos::RCP<Coupling::Adapter::MatrixRowColTransform> sggtransform_;
    //! transform object for structure interface matrix \f$S_{\Gamma I}\f$
    Teuchos::RCP<Coupling::Adapter::MatrixRowTransform> sgitransform_;
    //! transform object for structure interface matrix \f$S_{I \Gamma}\f$
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> sigtransform_;
    //! transform object for structure coupling matrix \f$C_{\Gamma \Gamma}^S\f$
    Teuchos::RCP<Coupling::Adapter::MatrixRowTransform> csggtransform_;
    //! transform object for fluid coupling matrix \f$C_{\Gamma \Gamma}^G\f$
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> cfggtransform_;
    //! transform object for structure coupling matrix \f$C_{\Gamma I}^S\f$
    Teuchos::RCP<Coupling::Adapter::MatrixRowTransform> csgitransform_;
    //! transform object for fluid coupling matrix \f$C_{I \Gamma}^F\f$
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> cfigtransform_;

    //!@}

    //! @name Some quantities to recover the Langrange multiplier at the end of each time step

    //! block \f$S_{\Gamma I,i+1}\f$ of structural matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const Core::LinAlg::SparseOperator> sgicur_;

    //! block \f$S_{\Gamma\Gamma,i+1}\f$ of structural matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const Core::LinAlg::SparseOperator> sggcur_;

    //! block \f$S_{\Gamma\Gamma,i+1}\f$ of structural matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const Core::LinAlg::SparseOperator> cgicur_;

    //! block \f$S_{\Gamma\Gamma,i+1}\f$ of structural matrix at current iteration \f$i+1\f$
    Teuchos::RCP<const Core::LinAlg::SparseOperator> cggcur_;

    //!@}
  };

}  // namespace PoroElast

FOUR_C_NAMESPACE_CLOSE

#endif
