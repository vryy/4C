/*----------------------------------------------------------------------*/
/*!
\file fsi_lungmonolithic_structuresplit.H
\brief Volume-coupled FSI (structure-split)

\maintainer Matthias Mayr

\level 3
*/
/*----------------------------------------------------------------------*/



#ifndef FSI_LUNG_MONOLITHICSTRUCTURESPLIT_H
#define FSI_LUNG_MONOLITHICSTRUCTURESPLIT_H

#include "fsi_lungmonolithic.H"
#include "fsi_matrixtransform.H"

namespace FSI
{
  /// monolithic FSI algorithm with overlapping interface equations
  /// for simulation of a specific class of bio problems (FSI airway
  /// model with attached balloon built of lung parenchyma)
  /*!
    \author lw
    \date 04/09
   */
  class LungMonolithicStructureSplit : public LungMonolithic
  {
   public:
    explicit LungMonolithicStructureSplit(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /*! do the setup for the monolithic system
    <pre>
    1.) setup coupling; right now, we use matching meshes at the interface
    2.) create combined map
    3.) create block system matrix
    </pre>
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

    UTILS::MatrixLogicalSplitAndTransform siitransform_;
    UTILS::MatrixLogicalSplitAndTransform sggtransform_;
    UTILS::MatrixLogicalSplitAndTransform sgitransform_;
    UTILS::MatrixLogicalSplitAndTransform sigtransform_;
    UTILS::MatrixColTransform aigtransform_;
    UTILS::MatrixColTransform aiGtransform_;

    UTILS::MatrixColTransform fmiitransform_;
    UTILS::MatrixColTransform fmGitransform_;
    UTILS::MatrixColTransform fmgitransform_;
    UTILS::MatrixColTransform fmiGtransform_;
    UTILS::MatrixColTransform fmGGtransform_;
    UTILS::MatrixColTransform fmgGtransform_;
    UTILS::MatrixColTransform addfmGGtransform_;

    /// split of constraint matrices
    UTILS::MatrixLogicalSplitAndTransform sciitransform_;
    UTILS::MatrixLogicalSplitAndTransform scgitransform_;
    UTILS::MatrixLogicalSplitAndTransform csiitransform_;
    UTILS::MatrixLogicalSplitAndTransform csigtransform_;
    UTILS::MatrixColTransform caiGtransform_;
  };
}  // namespace FSI

#endif
