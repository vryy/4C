/*----------------------------------------------------------------------*/
/*! \file

\brief State class for (in)stationary XFEM fluid problems involving embedded
fluid meshes

\level 2

 */
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_XFLUID_FLUID_STATE_HPP
#define FOUR_C_FLUID_XFLUID_FLUID_STATE_HPP

#include "4C_config.hpp"

#include "4C_fluid_xfluid_state.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Core::LinAlg
{
  class SparseOperator;
}

namespace XFEM
{
  class ConditionManager;
}

namespace FLD
{
  namespace UTILS
  {
    class XFluidFluidMapExtractor;
  }

  /**
   * Container class for the merged state vectors and maps of the intersected background
   * fluid and the embedded (ALE-)fluid.
   */
  class XFluidFluidState : public XFluidState
  {
   public:
    /// ctor
    explicit XFluidFluidState(const Teuchos::RCP<XFEM::ConditionManager>& condition_manager,
        const Teuchos::RCP<Core::Geo::CutWizard>& wizard,
        const Teuchos::RCP<XFEM::XFEMDofSet>& dofset,
        const Teuchos::RCP<const Epetra_Map>& xfluiddofrowmap,
        const Teuchos::RCP<const Epetra_Map>& xfluiddofcolmap,
        const Teuchos::RCP<const Epetra_Map>& embfluiddofrowmap);

    /// setup map extractors for dirichlet maps & velocity/pressure maps
    void SetupMapExtractors(const Teuchos::RCP<Discret::Discretization>& xfluiddiscret,
        const Teuchos::RCP<Discret::Discretization>& embfluiddiscret, const double& time);

    /// build merged fluid dirichlet map extractor
    void create_merged_dbc_map_extractor(
        Teuchos::RCP<const Core::LinAlg::MapExtractor> embfluiddbcmaps);

    //! @name Accessors
    //@{

    Teuchos::RCP<Core::LinAlg::MapExtractor> DBCMapExtractor() override { return xffluiddbcmaps_; }

    Teuchos::RCP<Core::LinAlg::MapExtractor> VelPresSplitter() override
    {
      return xffluidvelpressplitter_;
    }

    bool Destroy() override;

    Teuchos::RCP<Core::LinAlg::SparseMatrix> SystemMatrix() override;
    Teuchos::RCP<Epetra_Vector>& Residual() override { return xffluidresidual_; }
    Teuchos::RCP<Epetra_Vector>& Zeros() override { return xffluidzeros_; }
    Teuchos::RCP<Epetra_Vector>& IncVel() override { return xffluidincvel_; }
    Teuchos::RCP<Epetra_Vector>& Velnp() override { return xffluidvelnp_; }
    //@}

    void complete_coupling_matrices_and_rhs() override;

    //@name Map of the merged system
    //@{
    /// combined background and embedded fluid dof-map
    Teuchos::RCP<Epetra_Map> xffluiddofrowmap_;
    //@}

    //@name Map extractors of the merged system
    //@{
    /// extractor used for splitting fluid and embedded fluid
    Teuchos::RCP<FLD::UTILS::XFluidFluidMapExtractor> xffluidsplitter_;
    /// extractor used for splitting between velocity and pressure dof from the combined background
    /// & embedded fluid dof-map
    Teuchos::RCP<Core::LinAlg::MapExtractor> xffluidvelpressplitter_;
    /// combined background and embedded fluid map extractor for dirichlet-constrained dof
    Teuchos::RCP<Core::LinAlg::MapExtractor> xffluiddbcmaps_;
    //@}

    /// full system matrix for coupled background and embedded fluid
    Teuchos::RCP<Core::LinAlg::SparseOperator> xffluidsysmat_;

    /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> xffluidzeros_;

    /// (standard) residual vector (rhs for the incremental form),
    Teuchos::RCP<Epetra_Vector> xffluidresidual_;

    //! @name combined background and embedded fluid velocity and pressure at time n+1, n and
    //! increment
    //@{
    /// \f$ \mathbf{u}^{b\cup e,n+1} \f$
    Teuchos::RCP<Epetra_Vector> xffluidvelnp_;
    /// \f$ \mathbf{u}^{b\cup e,n+1} \f$
    Teuchos::RCP<Epetra_Vector> xffluidveln_;
    /// \f$ \Delta \mathbf{u}^{b\cup e,n+1}_{i+1} \f$
    Teuchos::RCP<Epetra_Vector> xffluidincvel_;
    //@}

   private:
    /// initialize all state members based on the merged fluid dof-rowmap
    void init_state_vectors();

    /// initialize the system matrix of the intersected fluid
    void init_system_matrix();

    /// embedded fluid dof-map
    Teuchos::RCP<const Epetra_Map> embfluiddofrowmap_;
  };

}  // namespace FLD

FOUR_C_NAMESPACE_CLOSE

#endif
