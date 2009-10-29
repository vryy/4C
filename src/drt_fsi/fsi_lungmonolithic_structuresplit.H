#ifdef CCADISCRET

#ifndef FSI_LUNG_MONOLITHICSTRUCTURESPLIT_H
#define FSI_LUNG_MONOLITHICSTRUCTURESPLIT_H

#include "fsi_lungmonolithic.H"
#include "fsi_monolithicstructuresplit.H"
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
  class LungMonolithicStructureSplit : public LungMonolithic
  {
  public:
    explicit LungMonolithicStructureSplit(Epetra_Comm& comm);

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

    UTILS::MatrixRowColTransform sggtransform_;
    UTILS::MatrixRowTransform sgitransform_;
    UTILS::MatrixColTransform sigtransform_;
    UTILS::MatrixColTransform aigtransform_;
    UTILS::MatrixColTransform aiGtransform_;

    UTILS::MatrixColTransform fmiitransform_;
    UTILS::MatrixColTransform fmGitransform_;
    UTILS::MatrixColTransform fmgitransform_;
    UTILS::MatrixColTransform fmiGtransform_;
    UTILS::MatrixColTransform fmGGtransform_;
    UTILS::MatrixColTransform fmgGtransform_;
    UTILS::MatrixColTransform addfmGGtransform_;
    UTILS::MatrixColTransform addfmgGtransform_;

    /// split of constraint matrices
    UTILS::MatrixRowTransform scgitransform_;
    UTILS::MatrixColTransform csigtransform_;
    UTILS::MatrixColTransform caiGtransform_;
  };
}

#endif
#endif



