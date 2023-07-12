/*----------------------------------------------------------------------*/
/*! \file

\brief State class for (in)stationary XFEM fluid problems involving embedded
fluid meshes

\level 2

 */
/*----------------------------------------------------------------------*/

#ifndef FLUID_XFLUID_FLUID_STATE_H
#define FLUID_XFLUID_FLUID_STATE_H

#include "fluid_xfluid_state.H"

// forward declarations
namespace CORE::LINALG
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
        const Teuchos::RCP<CORE::GEO::CutWizard>& wizard,
        const Teuchos::RCP<XFEM::XFEMDofSet>& dofset,
        const Teuchos::RCP<const Epetra_Map>& xfluiddofrowmap,
        const Teuchos::RCP<const Epetra_Map>& xfluiddofcolmap,
        const Teuchos::RCP<const Epetra_Map>& embfluiddofrowmap);

    /// dtor
    virtual ~XFluidFluidState() {}

    /// setup map extractors for dirichlet maps & velocity/pressure maps
    void SetupMapExtractors(const Teuchos::RCP<DRT::Discretization>& xfluiddiscret,
        const Teuchos::RCP<DRT::Discretization>& embfluiddiscret, const double& time);

    /// build merged fluid dirichlet map extractor
    void CreateMergedDBCMapExtractor(
        Teuchos::RCP<const CORE::LINALG::MapExtractor> embfluiddbcmaps);

    //! @name Accessors
    //@{

    Teuchos::RCP<CORE::LINALG::MapExtractor> DBCMapExtractor() { return xffluiddbcmaps_; }

    Teuchos::RCP<CORE::LINALG::MapExtractor> VelPresSplitter() { return xffluidvelpressplitter_; }

    virtual bool Destroy();

    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix();
    Teuchos::RCP<Epetra_Vector>& Residual() { return xffluidresidual_; }
    Teuchos::RCP<Epetra_Vector>& Zeros() { return xffluidzeros_; }
    Teuchos::RCP<Epetra_Vector>& IncVel() { return xffluidincvel_; }
    Teuchos::RCP<Epetra_Vector>& Velnp() { return xffluidvelnp_; }
    //@}

    void CompleteCouplingMatricesAndRhs();

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
    Teuchos::RCP<CORE::LINALG::MapExtractor> xffluidvelpressplitter_;
    /// combined background and embedded fluid map extractor for dirichlet-constrained dof
    Teuchos::RCP<CORE::LINALG::MapExtractor> xffluiddbcmaps_;
    //@}

    /// full system matrix for coupled background and embedded fluid
    Teuchos::RCP<CORE::LINALG::SparseOperator> xffluidsysmat_;

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
    void InitStateVectors();

    /// initialize the system matrix of the intersected fluid
    void InitSystemMatrix();

    /// embedded fluid dof-map
    Teuchos::RCP<const Epetra_Map> embfluiddofrowmap_;
  };

}  // namespace FLD

#endif
