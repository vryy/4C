/*----------------------------------------------------------------------*/
/*!
\file fsi_monolithicstructuresplit.H

\brief Solve FSI problem with matching grids using a monolithic scheme
with condensed structure interface displacements

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-15262
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef FSI_MONOLITHICSTRUCTURESPLIT_H
#define FSI_MONOLITHICSTRUCTURESPLIT_H

#include "fsi_monolithic.H"
#include "../drt_inpar/inpar_fsi.H"

// forward declarations
namespace ADAPTER
{
  class Structure;
  class Fluid;
}

namespace NOX
{
  namespace FSI
  {
    class AdaptiveNewtonNormF;
    class Group;
  }
}

namespace LINALG
{
  class BlockSparseMatrixBase;
}

namespace FSI
{
  class OverlappingBlockMatrix;

  namespace UTILS
  {
    class MatrixRowTransform;
    class MatrixColTransform;
    class MatrixRowColTransform;
  }
}

namespace FSI
{
  /// monolithic FSI algorithm with overlapping interface equations
  /*!

    Here the structural matrix is split whereas the fluid matrix is taken as
    it is.

    \sa MonolithicOverlap
    \author u.kue
    \date 06/08
   */
  class MonolithicStructureSplit : public BlockMonolithic
  {

    friend class FSI::FSIResultTest;

  public:
    explicit MonolithicStructureSplit(const Epetra_Comm& comm,
                                      const Teuchos::ParameterList& timeparams);

    /*! do the setup for the monolithic system

    <pre>

    1.) setup coupling; right now, we use matching meshes at the interface
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

    /// recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the end of each time step
    /// (i.e. condensed forces onto the structure) needed for rhs in next time step
    virtual void RecoverLagrangeMultiplier();

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig);

    /// the composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const;

    //! @name Methods for infnorm-scaling of the system

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b);

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b);

    //@}

    /// Output routine accounting for Lagrange multiplier at the interface
    virtual void Output();

    /// read restart data
    virtual void ReadRestart(int step);

    /// return Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface
    Teuchos::RCP<Epetra_Vector> GetLambda (){return lambda_; };

  protected:

    /// setup solver for global block system
    Teuchos::RCP<NOX::Epetra::LinearSystem>
    CreateLinearSystem(Teuchos::ParameterList& nlParams,
                       NOX::Epetra::Vector& noxSoln,
                       Teuchos::RCP<NOX::Utils> utils);

    /// setup of NOX convergence tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams,
                     Teuchos::RCP<NOX::Epetra::Group> grp);

    //! Extract the three field vectors from a given composed vector
    //!
    //! The condensed ale degrees of freedom have to be recovered
    //! from the fluid solution using a suitable velocity-displacement
    //! conversion. After that, the ale solution increment is projected onto the
    //! structure where a possible structural predictor has to
    //! be considered.
    //!
    //! We are dealing with NOX here, so we get absolute values. x is the sum of
    //! all increments up to this point.
    //!
    //! \sa  ADAPTER::FluidFSI::VelocityToDisplacement()
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
        const double fluidscale                 ///< residual scaling for fluid
    );

    /// block system matrix
    Teuchos::RCP<OverlappingBlockMatrix> systemmatrix_;

    /// coupling of fluid and ale at the free surface
    Teuchos::RCP<ADAPTER::Coupling> fscoupfa_;

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

    /// @name Some quantities to recover the Lagrange multiplier at the end of each time step

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the structure)
    //! evaluated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> lambda_;

    //! inner structural displacement increment \f$\Delta(\Delta d_{I,i+1}^{n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddiinc_;

    //! inner displacement solution of the structure at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> soliprev_;

    //! structural interface displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddginc_;

    //! fluid interface velocity increment \f$\Delta(\Delta u_{\Gamma,i+1}^{n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> duginc_;

    //! interface displacement solution of the structure at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> disgprev_;

    //! interface velocity solution of the fluid at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> velgprev_;

    //! block \f$S_{\Gamma I,i+1}\f$ of structural matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> sgicur_;

    //! block \f$S_{\Gamma I,i}\f$ of structural matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> sgiprev_;

    //! block \f$S_{\Gamma\Gamma,i+1}\f$ of structural matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> sggcur_;

    //! block \f$S_{\Gamma\Gamma,i}\f$ of structural matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> sggprev_;

    //@}

    /// additional ale residual to avoid incremental ale errors
    Teuchos::RCP<Epetra_Vector> aleresidual_;

    /// preconditioned block Krylov or block Gauss-Seidel linear solver
    INPAR::FSI::LinearBlockSolver linearsolverstrategy_;
  };
}

#endif
