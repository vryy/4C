#ifdef CCADISCRET

#ifndef FSI_MONOLITHICOVERLAP_H
#define FSI_MONOLITHICOVERLAP_H

#include <NOX_Direction_UserDefinedFactory.H>

#include "fsi_monolithic.H"
#include "fsi_overlapprec.H"
#include "fsi_matrixtransform.H"

namespace NOX
{
namespace FSI
{
  class PartialNormF;
}
}

namespace FSI
{

  /// monolithic FSI algorithm with overlapping interface equations
  /*!

    This is research. Highly experimental. Subject to change.

    Combine structure, fluid and ale field in one huge block matrix. Matching
    nodes. Overlapping equations at the interface.

    The structure equations contain both internal and interface part. Fluid
    and ale blocks are reduced to their respective internal parts. The fluid
    interface equations are added to the structure interface equations. There
    are no ale equations at the interface.

    Based on Newton's method within NOX. NOX computes the sum of all Newton
    increments. The evaluation method computeF() is always called with
    the sum x. However the meaning of this sum depends on the field blocks
    used.

    - The structure block calculates the displacement increments:
      \f$ \Delta \mathbf{d}^{n+1}_{i+1} = \mathbf{d}^{n+1}_{i+1} - \mathbf{d}^{n} \f$

    - The fluid block calculates the velocity (and pressure) increments:
      \f$ \Delta \mathbf{u}^{n+1}_{i+1} = \mathbf{u}^{n+1}_{i+1} - \mathbf{u}^{n} \f$

    - The ale block calculates the absolute mesh displacement:
      \f$ \Delta \mathbf{d}^{G,n+1} = \Delta \mathbf{d}^{n+1}_{i+1} + \mathbf{d}^{n} \f$

    We assume a very simple velocity -- displacement relation at the interface
    \f$ \mathbf{u}^{n+1}_{i+1} = \frac{1}{\Delta t} \Delta \mathbf{d}^{n+1}_{i+1} \f$

    \author u.kue
    \date 02/08
   */
  class MonolithicOverlap : public Monolithic,
                            public NOX::Direction::UserDefinedFactory
  {
  public:
    explicit MonolithicOverlap(Epetra_Comm& comm);

    //! @name NOX methods

    /// compute FSI residual
    bool computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag);

    /// compute FSI block matrix
    bool computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac);

    /// preconditioner
    bool computePreconditioner(const Epetra_Vector &x, Epetra_Operator &M, Teuchos::ParameterList *precParams=0);

    //@}

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Epetra_Vector& f);

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat);

    //@}

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    /// the composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const { return systemmatrix_; }

    /// create my own direction object
    /*!
      MonolithicOverlap is a (inherits from)
      NOX::Direction::UserDefinedFactory. This is an implementation
      detail. This way we can construct a specialized direction object at a
      place where we know about the status tests. This is the whole point
      here. Our specialized direction is of the type NOX::FSI::Newton, the
      normal Newton direction enhanced with adaptive tolerance control for the
      internal linear (iterative) solver.
     */
    virtual Teuchos::RCP<NOX::Direction::Generic> buildDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
                                                                 Teuchos::ParameterList& params) const;

  protected:

    /// setup solver for global block system
    Teuchos::RCP<NOX::Epetra::LinearSystem>
    CreateLinearSystem(ParameterList& nlParams,
                       NOX::Epetra::Vector& noxSoln,
                       Teuchos::RCP<NOX::Utils> utils);

    /// setup of NOX convergance tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Epetra::Group> grp);

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


  private:

    /// setup list with default parameters
    void SetDefaultParameters(const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list);

    /// build block vector from field vectors
    void SetupVector(Epetra_Vector &f,
                     Teuchos::RCP<const Epetra_Vector> sv,
                     Teuchos::RCP<const Epetra_Vector> fv,
                     Teuchos::RCP<const Epetra_Vector> av,
                     double fluidscale);

    ///
    void AddFluidInterface(double scale,
                           const LINALG::SparseMatrix& fgg,
                           LINALG::SparseMatrix& s);

    /// block system matrix
    Teuchos::RCP<OverlappingBlockMatrix> systemmatrix_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    UTILS::MatrixRowColTransform fggtransform_;
    UTILS::MatrixRowTransform fgitransform_;
    UTILS::MatrixColTransform figtransform_;
    UTILS::MatrixColTransform aigtransform_;
    UTILS::MatrixColTransform fmiitransform_;
    UTILS::MatrixRowColTransform fmgitransform_;

    ///@}

    /// keep the status tests available so we can connect them with our
    /// adaptive Newton direction
    std::vector<Teuchos::RCP<NOX::FSI::PartialNormF> > statustests_;
  };
}

#endif
#endif
