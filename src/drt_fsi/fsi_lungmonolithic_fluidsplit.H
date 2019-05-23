/*----------------------------------------------------------------------*/
/*!
\brief Volume-coupled FSI (fluid-split)

\maintainer Matthias Mayr

\level 3
*/
/*----------------------------------------------------------------------*/



#ifndef FSI_LUNG_MONOLITHIC_FLUIDSPLIT_H
#define FSI_LUNG_MONOLITHIC_FLUIDSPLIT_H

#include "fsi_lungmonolithic.H"

namespace FSI
{
  // forward declarations
  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }  // namespace UTILS

  /// monolithic FSI algorithm with overlapping interface equations
  /// for simulation of a specific class of bio problems (FSI airway
  /// model with attached balloon built of lung parenchyma)
  /*!
    \author lw
    \date 04/09
   */
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
    virtual void SetupSystem();

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat);

    //@}

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

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
    virtual void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
        Teuchos::RCP<const Epetra_Vector>& ax);

    /// build block vector from field vectors
    virtual void SetupVector(Epetra_Vector& f,
        Teuchos::RCP<const Epetra_Vector> sv,  ///< structure vector
        Teuchos::RCP<const Epetra_Vector> fv,  ///< fluid vector
        Teuchos::RCP<const Epetra_Vector> av,  ///< ale vector
        Teuchos::RCP<const Epetra_Vector> cv,  ///< constraint vector
        double fluidscale);                    ///< scaling

   private:
    /*! \brief Create the combined DOF row map for the FSI problem
     *
     *  Combine the DOF row maps of structure, fluid and ALE to an global FSI
     *  DOF row map.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void CreateCombinedDofRowMap();

    /*! \brief Setup the Dirichlet map extractor
     *
     *  Create a map extractor #dbcmaps_ for the Dirichlet degrees of freedom
     *  for the entire FSI problem. This is done just by combining the
     *  condition maps and other maps from structure, fluid and ALE to a FSI-global
     *  condition map and other map.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void SetupDBCMapExtractor() { dserror("Not implemented, yet."); }

    /// setup RHS contributions based on single field residuals
    void SetupRHSResidual(Epetra_Vector& f);

    /// setup RHS contributions based on the Lagrange multiplier field
    void SetupRHSLambda(Epetra_Vector& f);

    /// setup RHS contributions based on terms for first nonlinear iteration
    void SetupRHSFirstiter(Epetra_Vector& f);

    /// transformation of fluid matrix
    Teuchos::RCP<UTILS::MatrixRowColTransform> fggtransform_;
    Teuchos::RCP<UTILS::MatrixRowTransform> fgitransform_;
    Teuchos::RCP<UTILS::MatrixRowTransform> fgGtransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> figtransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> fGgtransform_;

    /// transformation of shape derivative matrix
    Teuchos::RCP<UTILS::MatrixColTransform> fmiitransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> fmGitransform_;
    Teuchos::RCP<UTILS::MatrixRowColTransform> fmgitransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> fmigtransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> fmGgtransform_;
    Teuchos::RCP<UTILS::MatrixRowColTransform> fmggtransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> fmiGtransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> fmGGtransform_;
    Teuchos::RCP<UTILS::MatrixRowColTransform> fmgGtransform_;

    /// transformation of additional shape derivative matrix (volume constraint)
    Teuchos::RCP<UTILS::MatrixColTransform> addfmGGtransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> addfmGgtransform_;

    /// transformation of fluid constraint matrix
    Teuchos::RCP<UTILS::MatrixRowTransform> fcgitransform_;

    /// transformation of ale matrix
    Teuchos::RCP<UTILS::MatrixColTransform> aigtransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> aiGtransform_;

    /// transformation of constraint "ale" matrix
    Teuchos::RCP<UTILS::MatrixColTransform> caiGtransform_;
  };
}  // namespace FSI

#endif
