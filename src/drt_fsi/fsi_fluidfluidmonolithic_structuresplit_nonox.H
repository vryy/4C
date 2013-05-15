

#ifndef FSI_FLUIDFLUIDMONOLITHICSTRUCTURESPLIT_H
#define FSI_FLUIDFLUIDMONOLITHICSTRUCTURESPLIT_H

#include "../drt_inpar/inpar_xfem.H"
#include "fsi_monolithic_nonox.H"

// forward declarations
namespace ADAPTER
{
  class Coupling;
}

namespace FSI
{
  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }
}

namespace FSI
{

  /// monolithic Fluid-Fluid FSI algorithm Structure Split
  /*!

    Here the structural matrix is split whereas the fluid matrix is taken as
    it is.

    \author Shadan Shahmiri
    \date  11/2011
  */
  class FluidFluidMonolithicStructureSplitNoNOX : public MonolithicNoNOX
  {
    friend class FSI::FSIResultTest;

  public:
    explicit FluidFluidMonolithicStructureSplitNoNOX(const Epetra_Comm& comm,
                                                const Teuchos::ParameterList& timeparams);

    /*! do the setup for the monolithic system

    <pre>

    1.) setup coupling; right now, we use matching meshes at the interface
    2.) create combined map
    3.) create block system matrix

    </pre>

    */
    virtual void SetupSystem();

    void SetupNewSystem();

    //! @name Apply current field state to system

    /// setup composed right hand side from field solvers
    virtual void SetupRHS(Epetra_Vector& f, bool firstcall);

    /// setup composed system matrix from field solvers
    //void SetupSystemMatrixold(LINALG::BlockSparseMatrixBase& mat);

    virtual void SetupSystemMatrix();

    //@}

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b);

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b);

    virtual Teuchos::RCP<Epetra_Map> CombinedDBCMap();

    // newton for xfem
    virtual void Newton();

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

    virtual void Update();

    virtual void ReadRestart(int step);

    virtual void Output();

    virtual void PrepareTimeStep();

    virtual void RecoverLagrangeMultiplier();

    virtual void BuildCovergenceNorms();


  private:

    /// build block vector from field vectors
    void SetupVector(Epetra_Vector &f,
                     Teuchos::RCP<const Epetra_Vector> sv,
                     Teuchos::RCP<const Epetra_Vector> fv,
                     Teuchos::RCP<const Epetra_Vector> av,
                     double fluidscale);

    /// block system matrix
    //Teuchos::RCP<OverlappingBlockMatrix> systemmatrix_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    Teuchos::RCP<UTILS::MatrixRowColTransform> sggtransform_;
    Teuchos::RCP<UTILS::MatrixRowTransform> sgitransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> sigtransform_;
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

    Teuchos::RCP<FLD::UTILS::FluidXFluidMapExtractor> xfluidfluidsplitter_;

    enum INPAR::XFEM::Monolithic_xffsi_Approach monolithic_approach_;

    int currentstep_;

    int relaxing_ale_every_;

    bool relaxing_ale_;

    /// @name Some quantities to recover the Langrange multiplier at the end of each time step

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the structure)
    //! evalutated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    // lambda lives at the slave side (here at stucture)
    Teuchos::RCP<Epetra_Vector> lambda_;

    //! interface force \f$f_{\Gamma,i+1}^{S,n+1}\f$ onto the structure at current NOX iteration \f$i+1\f$
    //Teuchos::RCP<const Epetra_Vector> fgcur_;

    //! interface force \f$f_{\Gamma,i}^{S,n+1}\f$ onto the structure at previous NOX iteration \f$i\f$
    //Teuchos::RCP<const Epetra_Vector> fgpre_;

    //! inner structural displacement increment \f$\Delta(\Delta d_{I,i+1}^{n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddiinc_;

    //! inner displacement solution of the structure at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> solipre_;

    //! structural interface displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddginc_;

    //! interface displacement solution of the structure at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> solgpre_;

    //! block \f$S_{\Gamma I,i+1}\f$ of structural matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> sgicur_;

    //! block \f$S_{\Gamma I,i}\f$ of structural matrix at previous NOX iteration \f$i\f$
    //Teuchos::RCP<const LINALG::SparseMatrix> sgipre_;

    //! block \f$S_{\Gamma\Gamma,i+1}\f$ of structural matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> sggcur_;

    //! block \f$S_{\Gamma\Gamma,i}\f$ of structural matrix at previous NOX iteration \f$i\f$
    //Teuchos::RCP<const LINALG::SparseMatrix> sggpre_;

    //! fluid velocity increment of interface DOFs due to predictor or inhomogeneous DBCs
    Teuchos::RCP<Epetra_Vector> dugpre_;


  };
}

#endif
