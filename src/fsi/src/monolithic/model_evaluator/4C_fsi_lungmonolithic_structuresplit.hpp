/*----------------------------------------------------------------------*/
/*! \file
\brief Volume-coupled FSI (structure-split)


\level 3
*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_FSI_LUNGMONOLITHIC_STRUCTURESPLIT_HPP
#define FOUR_C_FSI_LUNGMONOLITHIC_STRUCTURESPLIT_HPP

#include "4C_config.hpp"

#include "4C_fsi_lungmonolithic.hpp"
#include "4C_linalg_matrixtransform.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FSI
{
  /// monolithic FSI algorithm with overlapping interface equations
  /// for simulation of a specific class of bio problems (FSI airway
  /// model with attached balloon built of lung parenchyma)
  class LungMonolithicStructureSplit : public LungMonolithic
  {
   public:
    explicit LungMonolithicStructureSplit(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /*! do the setup for the monolithic system
      1.) setup coupling; right now, we use matching meshes at the interface
      2.) create combined map
      3.) create block system matrix
    */
    void SetupSystem() override;

    /// setup composed system matrix from field solvers
    void setup_system_matrix(Core::LinAlg::BlockSparseMatrixBase& mat) override;

    //@}

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

    Core::LinAlg::MatrixLogicalSplitAndTransform siitransform_;
    Core::LinAlg::MatrixLogicalSplitAndTransform sggtransform_;
    Core::LinAlg::MatrixLogicalSplitAndTransform sgitransform_;
    Core::LinAlg::MatrixLogicalSplitAndTransform sigtransform_;
    Core::LinAlg::MatrixColTransform aigtransform_;
    Core::LinAlg::MatrixColTransform ai_gtransform_;

    Core::LinAlg::MatrixColTransform fmiitransform_;
    Core::LinAlg::MatrixColTransform fm_gitransform_;
    Core::LinAlg::MatrixColTransform fmgitransform_;
    Core::LinAlg::MatrixColTransform fmi_gtransform_;
    Core::LinAlg::MatrixColTransform fm_g_gtransform_;
    Core::LinAlg::MatrixColTransform fmg_gtransform_;
    Core::LinAlg::MatrixColTransform addfm_g_gtransform_;

    /// split of constraint matrices
    Core::LinAlg::MatrixLogicalSplitAndTransform sciitransform_;
    Core::LinAlg::MatrixLogicalSplitAndTransform scgitransform_;
    Core::LinAlg::MatrixLogicalSplitAndTransform csiitransform_;
    Core::LinAlg::MatrixLogicalSplitAndTransform csigtransform_;
    Core::LinAlg::MatrixColTransform cai_gtransform_;
  };
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
