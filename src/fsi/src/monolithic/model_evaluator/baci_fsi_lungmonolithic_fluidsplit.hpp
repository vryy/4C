/*----------------------------------------------------------------------*/
/*! \file
\brief Volume-coupled FSI (fluid-split)


\level 3
*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_FSI_LUNGMONOLITHIC_FLUIDSPLIT_HPP
#define FOUR_C_FSI_LUNGMONOLITHIC_FLUIDSPLIT_HPP

#include "baci_config.hpp"

#include "baci_fsi_lungmonolithic.hpp"

BACI_NAMESPACE_OPEN

// forward declarations
namespace CORE::LINALG
{
  class MatrixRowTransform;
  class MatrixColTransform;
  class MatrixRowColTransform;
}  // namespace CORE::LINALG

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
    void SetupSystem() override;

    /// setup composed system matrix from field solvers
    void SetupSystemMatrix(CORE::LINALG::BlockSparseMatrixBase& mat) override;

    /// Extract initial guess from fields
    void InitialGuess(Teuchos::RCP<Epetra_Vector> ig) override;

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
    void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
        Teuchos::RCP<const Epetra_Vector>& ax) override;

    /// build block vector from field vectors
    void SetupVector(Epetra_Vector& f,
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
    void CreateCombinedDofRowMap() override;

    /*! \brief Setup the Dirichlet map extractor
     *
     *  Create a map extractor #dbcmaps_ for the Dirichlet degrees of freedom
     *  for the entire FSI problem. This is done just by combining the
     *  condition maps and other maps from structure, fluid and ALE to a FSI-global
     *  condition map and other map.
     */
    void SetupDBCMapExtractor() override { dserror("Not implemented, yet."); }

    /// setup RHS contributions based on single field residuals
    void SetupRHSResidual(Epetra_Vector& f) override;

    /// setup RHS contributions based on the Lagrange multiplier field
    void SetupRHSLambda(Epetra_Vector& f) override;

    /// setup RHS contributions based on terms for first nonlinear iteration
    void SetupRHSFirstiter(Epetra_Vector& f) override;

    /// transformation of fluid matrix
    Teuchos::RCP<CORE::LINALG::MatrixRowColTransform> fggtransform_;
    Teuchos::RCP<CORE::LINALG::MatrixRowTransform> fgitransform_;
    Teuchos::RCP<CORE::LINALG::MatrixRowTransform> fgGtransform_;
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> figtransform_;
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> fGgtransform_;

    /// transformation of shape derivative matrix
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> fmiitransform_;
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> fmGitransform_;
    Teuchos::RCP<CORE::LINALG::MatrixRowColTransform> fmgitransform_;
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> fmigtransform_;
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> fmGgtransform_;
    Teuchos::RCP<CORE::LINALG::MatrixRowColTransform> fmggtransform_;
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> fmiGtransform_;
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> fmGGtransform_;
    Teuchos::RCP<CORE::LINALG::MatrixRowColTransform> fmgGtransform_;

    /// transformation of additional shape derivative matrix (volume constraint)
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> addfmGGtransform_;
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> addfmGgtransform_;

    /// transformation of fluid constraint matrix
    Teuchos::RCP<CORE::LINALG::MatrixRowTransform> fcgitransform_;

    /// transformation of ale matrix
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> aigtransform_;
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> aiGtransform_;

    /// transformation of constraint "ale" matrix
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> caiGtransform_;
  };
}  // namespace FSI

BACI_NAMESPACE_CLOSE

#endif
