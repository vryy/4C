

#ifndef FSI_MORTARMONOLITHICSTRUCTURESPLIT_H
#define FSI_MORTARMONOLITHICSTRUCTURESPLIT_H

#include "../drt_inpar/inpar_fsi.H"
#include "fsi_monolithic.H"
#include "fsi_overlapprec.H"

// forward declarations
namespace ADAPTER
{
  class Coupling;
  class CouplingMortar;
}

namespace FSI
{
  namespace UTILS
  {
    class MatrixColTransform;
    class SlideAleUtils;
  }
}

namespace FSI
{

  /// monolithic FSI algorithm with overlapping non-matching interface equations
  /*!

    Here the structural matrix is split whereas the fluid matrix is taken as
    it is.

    \sa MortarMonolithicStructureSplit
    \author tk
    \date 06/08
   */
  class MortarMonolithicStructureSplit : public BlockMonolithic
  {
  public:
    explicit MortarMonolithicStructureSplit(const Epetra_Comm& comm,
                                            const Teuchos::ParameterList& timeparams);

    /*! do the setup for the monolithic system

    <pre>

    1.) setup coupling
    2.) create combined map
    3.) create block system matrix

    </pre>

    */
    virtual void SetupSystem();

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Epetra_Vector& f, bool firstcall=false);

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat);

    //@}

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    /// the composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const { return systemmatrix_; }

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b);

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b);

    /// read restart
    virtual void ReadRestart(int step);

    /// recover Lagrange multiplier at the interface (i.e. condensed forces onto the structure)
    /// needed for rhs in next time step
    virtual void RecoverLagrangeMultiplier(){};

  protected:

    virtual void Update();

    virtual void Output();

    /// setup solver for global block system
    Teuchos::RCP<NOX::Epetra::LinearSystem>
    CreateLinearSystem(Teuchos::ParameterList& nlParams,
                       NOX::Epetra::Vector& noxSoln,
                       Teuchos::RCP<NOX::Utils> utils);

    /// setup of NOX convergence tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Epetra::Group> grp);

    /// extract the three field vectors from a given composed vector
    /*!
      We are dealing with NOX here, so we get absolute values. x is the sum of
      all increments up to this point.
     */
    virtual void ExtractFieldVectors(
        Teuchos::RCP<const Epetra_Vector> x,    ///< composed vector that contains all field vectors
        Teuchos::RCP<const Epetra_Vector>& sx,  ///< structural displacements
        Teuchos::RCP<const Epetra_Vector>& fx,  ///< fluid velocities and pressure
        Teuchos::RCP<const Epetra_Vector>& ax   ///< ale displacements
    );


  private:

    /// build block vector from field vectors
    void SetupVector(
        Epetra_Vector &f,                       ///< composed vector that contains all field vectors
        Teuchos::RCP<const Epetra_Vector> sv,   ///< structural displacements
        Teuchos::RCP<const Epetra_Vector> fv,   ///< fluid velocities and pressure
        Teuchos::RCP<const Epetra_Vector> av,   ///< ale displacements
        double fluidscale                       ///< residual scaling for fluid
    );

    /// block system matrix
    Teuchos::RCP<OverlappingBlockMatrix> systemmatrix_;

    /// coupling of fluid and ale (interface only)
    Teuchos::RCP< ::ADAPTER::Coupling> icoupfa_;

    /// coupling of fluid and ale at the free surface
    Teuchos::RCP< ::ADAPTER::Coupling> fscoupfa_;

    /// coupling of structure and fluid at the interface
    Teuchos::RCP< ::ADAPTER::CouplingMortar>               coupsfm_;

    ///communicator
    const Epetra_Comm& comm_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    Teuchos::RCP<UTILS::MatrixColTransform> aigtransform_;

    Teuchos::RCP<UTILS::MatrixColTransform> fmiitransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> fmgitransform_;

    Teuchos::RCP<UTILS::MatrixColTransform> fsaigtransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> fsmgitransform_;

    ///@}

    /// @name infnorm scaling

    Teuchos::RCP<Epetra_Vector> srowsum_;
    Teuchos::RCP<Epetra_Vector> scolsum_;
    Teuchos::RCP<Epetra_Vector> arowsum_;
    Teuchos::RCP<Epetra_Vector> acolsum_;

    //@}

    /// additional ale residual to avoid incremental ale errors
    Teuchos::RCP<Epetra_Vector> aleresidual_;

    /// preconditioned block Krylov or block Gauss-Seidel linear solver
    INPAR::FSI::LinearBlockSolver linearsolverstrategy_;

    /// ale movement relative to structure (none, slide_curr, slide_ref)
    INPAR::FSI::SlideALEProj aleproj_;

    bool notsetup_;  ///< indicates if Setup has not been called yet

    Teuchos::RCP<FSI::UTILS::SlideAleUtils> slideale_;  ///< Sliding Ale helper class

    Teuchos::RCP<Epetra_Vector> iprojdispinc_;        ///<displacement of fluid side of the interface
    Teuchos::RCP<Epetra_Vector> iprojdisp_;        ///<displacement of fluid side of the interface

  };
}

#endif
