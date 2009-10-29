#ifdef CCADISCRET

#ifndef FSI_LUNG_MONOLITHIC_FLUIDSPLIT_H
#define FSI_LUNG_MONOLITHIC_FLUIDSPLIT_H

#include "fsi_lungmonolithic.H"
#include "fsi_matrixtransform.H"
#include "../drt_constraint/constraintdofset.H"

namespace FSI
{

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
    explicit LungMonolithicFluidSplit(Epetra_Comm& comm);

    /*! do the setup for the monolithic system
    <pre>
    1.) setup coupling; right now, we use matching meshes at the interface
    2.) create combined map
    3.) create block system matrix
    </pre>
    */
    virtual void SetupSystem();

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Epetra_Vector& f, bool firstcall=false);

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
                                     Teuchos::RCP<const Epetra_Vector>& sx,
                                     Teuchos::RCP<const Epetra_Vector>& fx,
                                     Teuchos::RCP<const Epetra_Vector>& ax);

    /// build block vector from field vectors
    virtual void SetupVector(Epetra_Vector &f,
                             Teuchos::RCP<const Epetra_Vector> sv,         ///< structure vector
                             Teuchos::RCP<const Epetra_Vector> fv,         ///< fluid vector
                             Teuchos::RCP<const Epetra_Vector> av,         ///< ale vector
                             Teuchos::RCP<const Epetra_Vector> cv,         ///< constraint vector
                             double fluidscale);                           ///< scaling


  private:

    /// transformation of fluid matrix
    UTILS::MatrixRowColTransform fggtransform_;
    UTILS::MatrixRowTransform fgitransform_;
    UTILS::MatrixRowTransform fgGtransform_;
    UTILS::MatrixColTransform figtransform_;
    UTILS::MatrixColTransform fGgtransform_;

    /// transformation of shape derivative matrix
    UTILS::MatrixColTransform fmiitransform_;
    UTILS::MatrixColTransform fmGitransform_;
    UTILS::MatrixRowColTransform fmgitransform_;
    UTILS::MatrixColTransform fmigtransform_;
    UTILS::MatrixColTransform fmGgtransform_;
    UTILS::MatrixRowColTransform fmggtransform_;
    UTILS::MatrixColTransform fmiGtransform_;
    UTILS::MatrixColTransform fmGGtransform_;
    UTILS::MatrixRowColTransform fmgGtransform_;

    /// transformation of additional shape derivative matrix (volume constraint)
    UTILS::MatrixRowColTransform addfmgGtransform_;
    UTILS::MatrixColTransform addfmGGtransform_;
    UTILS::MatrixRowColTransform addfmggtransform_;
    UTILS::MatrixColTransform addfmGgtransform_;

    /// transformation of fluid constraint matrix
    UTILS::MatrixRowTransform fcgitransform_;

    /// transformation of ale matrix
    UTILS::MatrixColTransform aigtransform_;
    UTILS::MatrixColTransform aiGtransform_;

    /// transformation of constraint fluid matrix
    UTILS::MatrixColTransform cfigtransform_;

    /// transformation of constraint "ale" matrix
    UTILS::MatrixColTransform caigtransform_;
    UTILS::MatrixColTransform caiGtransform_;

  };
}

#endif
#endif



