/*----------------------------------------------------------------------*/
/*! \file
\brief Volume-coupled FSI (fluid-split)


\level 3
*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_FSI_LUNGMONOLITHIC_FLUIDSPLIT_HPP
#define FOUR_C_FSI_LUNGMONOLITHIC_FLUIDSPLIT_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter_converter.hpp"
#include "4C_fsi_lungmonolithic.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FSI
{
  /// monolithic FSI algorithm with overlapping interface equations
  /// for simulation of a specific class of bio problems (FSI airway
  /// model with attached balloon built of lung parenchyma)
  class LungMonolithicFluidSplit : public LungMonolithic
  {
   public:
    explicit LungMonolithicFluidSplit(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /*! do the setup for the monolithic system
      1.) setup coupling; right now, we use matching meshes at the interface
      2.) create combined map
      3.) create block system matrix
    */
    void setup_system() override;

    /// setup composed system matrix from field solvers
    void setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat) override;

    /// Extract initial guess from fields
    void initial_guess(Teuchos::RCP<Epetra_Vector> ig) override;

   protected:
    /// extract the three field vectors from a given composed vector
    /*!
      We are dealing with NOX here, so we get absolute values. x is the sum of
      all increments up to this point.

      \param x  (i) composed vector that contains all field vectors
      \param sx (o) structural displacements
      \param fx (o) fluid velocities and pressure
      \param ax (o) ale displacements
    */
    void extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
        Teuchos::RCP<const Epetra_Vector>& ax) override;

    /// build block vector from field vectors
    void setup_vector(Epetra_Vector& f,
        Teuchos::RCP<const Epetra_Vector> sv,  ///< structure vector
        Teuchos::RCP<const Epetra_Vector> fv,  ///< fluid vector
        Teuchos::RCP<const Epetra_Vector> av,  ///< ale vector
        Teuchos::RCP<const Epetra_Vector> cv,  ///< constraint vector
        double fluidscale) override;           ///< scaling

   private:
    /*! \brief Create the combined DOF row map for the FSI problem
     *
     *  Combine the DOF row maps of structure, fluid and ALE to an global FSI
     *  DOF row map.
     */
    void create_combined_dof_row_map() override;

    /*! \brief Setup the Dirichlet map extractor
     *
     *  Create a map extractor #dbcmaps_ for the Dirichlet degrees of freedom
     *  for the entire FSI problem. This is done just by combining the
     *  condition maps and other maps from structure, fluid and ALE to a FSI-global
     *  condition map and other map.
     */
    void setup_dbc_map_extractor() override { FOUR_C_THROW("Not implemented, yet."); }

    /// setup RHS contributions based on single field residuals
    void setup_rhs_residual(Epetra_Vector& f) override;

    /// setup RHS contributions based on the Lagrange multiplier field
    void setup_rhs_lambda(Epetra_Vector& f) override;

    /// setup RHS contributions based on terms for first nonlinear iteration
    void setup_rhs_firstiter(Epetra_Vector& f) override;

    /// transformation of fluid matrix
    Teuchos::RCP<Coupling::Adapter::MatrixRowColTransform> fggtransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixRowTransform> fgitransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixRowTransform> fg_gtransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> figtransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> f_ggtransform_;

    /// transformation of shape derivative matrix
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> fmiitransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> fm_gitransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixRowColTransform> fmgitransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> fmigtransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> fm_ggtransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixRowColTransform> fmggtransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> fmi_gtransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> fm_g_gtransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixRowColTransform> fmg_gtransform_;

    /// transformation of additional shape derivative matrix (volume constraint)
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> addfm_g_gtransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> addfm_ggtransform_;

    /// transformation of fluid constraint matrix
    Teuchos::RCP<Coupling::Adapter::MatrixRowTransform> fcgitransform_;

    /// transformation of ale matrix
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> aigtransform_;
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> ai_gtransform_;

    /// transformation of constraint "ale" matrix
    Teuchos::RCP<Coupling::Adapter::MatrixColTransform> cai_gtransform_;
  };
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
