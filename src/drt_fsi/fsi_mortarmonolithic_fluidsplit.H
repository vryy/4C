

#ifndef FSI_MORTARMONOLITHICFLUID_H
#define FSI_MORTARMONOLITHICFLUID_H

#include "fsi_monolithic.H"
#include "fsi_overlapprec.H"
#include "../drt_inpar/inpar_fsi.H"

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

  Here the fluid matrix is split whereas the structure matrix is taken as
  it is.

  \sa MortarMonolithicFluidSplit
  \author tk
  \date 06/08
   */
  class MortarMonolithicFluidSplit : public BlockMonolithic
  {
  public:
    explicit MortarMonolithicFluidSplit(const Epetra_Comm& comm,
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

    /// start a new time step
    virtual void PrepareTimeStep();

    /// recover Lagrange multiplier at the interface needed for rhs in next time step
    virtual void RecoverLagrangeMultiplier();

  protected:

    virtual void Update();

    virtual void Output();

    /// setup solver for global block system
    Teuchos::RCP<NOX::Epetra::LinearSystem>
    CreateLinearSystem(ParameterList& nlParams,
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
    /* Here, all contributions to rhs vector, that are used in every Newton iteration, are assembled.
     *  a) residuals from single fields
     *  b) Lagrange multiplier \f$\lambda^n\f$ from last time step
     */
    void SetupVector(
        Epetra_Vector &f,                       ///< composed vector that contains all field vectors
        Teuchos::RCP<const Epetra_Vector> sv,   ///< structural displacements
        Teuchos::RCP<const Epetra_Vector> fv,   ///< fluid velocities and pressure
        Teuchos::RCP<const Epetra_Vector> av,   ///< ale displacements
        double fluidscale                       ///< residual scaling for fluid
    );

    /// block system matrix
    /* System matrix has a 4x4-block structure corresponding to the vector of unknowns
     *
     * \f$\Delta x^T = [\Delta d_I^{S,n+1}~\Delta d_\Gamma^{S,n+1}~\Delta u_I^{F,n+1}~\Delta d_I^{G,n+1}]\f$
     */
    Teuchos::RCP<OverlappingBlockMatrix> systemmatrix_;

    ///communicator
    const Epetra_Comm& comm_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    /// coupling of fluid and ale (interface only)
    Teuchos::RCP< ::ADAPTER::Coupling> icoupfa_;

    /// coupling of fluid and ale at the free surface
    Teuchos::RCP< ::ADAPTER::Coupling> fscoupfa_;

    /// coupling of structure and fluid at the interface
    Teuchos::RCP< ::ADAPTER::CouplingMortar>               coupsfm_;

    Teuchos::RCP<UTILS::MatrixColTransform> aigtransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> fmiitransform_;

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

    /// @name Recovery of Langrange multiplier at the end of each time step

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the fluid)
    //! evalutated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> lambda_;

    //! interface force \f$f_{\Gamma,i+1}^{F,n+1}\f$ onto the fluid at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const Epetra_Vector> fgcur_;

    //! interface force \f$f_{\Gamma,i}^{F,n+1}\f$ onto the fluid at previous NOX iteration \f$i\f$
    Teuchos::RCP<const Epetra_Vector> fgpre_;

    //! interface structure displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddginc_;

    //! inner fluid velocity increment \f$\Delta(\Delta u_{I,i+1}^{n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> duiinc_;

    //! interface displacement solution of the structure at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> disgpre_;

    //! inner velocity solution of fluid at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> velipre_;

    //! interface velocity solution of the fluid at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> velgpre_;

    //! inner ALE displacement solution at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> aleipre_;

    //! interface ALE displacement solution at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> alegpre_;

    //! inner ALE displacement increment \f$\Delta(\Delta d_{I,i+1}^{G,n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddialeinc_;

    //! block \f$F_{\Gamma I,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fgicur_;

    //! block \f$F_{\Gamma I,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fgipre_;

    //! block \f$F_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fggcur_;

    //! block \f$F_{\Gamma\Gamma,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fggpre_;

    //! block \f$F_{\Gamma I,i+1}^G\f$ of ALE matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmgicur_;

    //! block \f$F_{\Gamma I,i}^G\f$ of ALE matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmgipre_;

    //! block \f$F_{\Gamma\Gamma,i+1}^G\f$ of ALE matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmggcur_;

    //! block \f$F_{\Gamma\Gamma,i}^G\f$ of ALE matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmggpre_;

    //@}

  };
}

#endif
