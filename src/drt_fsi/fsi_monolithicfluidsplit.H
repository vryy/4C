/*----------------------------------------------------------------------*/
/*!
\file fsi_monolithicfluidsplit.H

\brief Solve FSI problem with matching grids using a monolithic scheme
with condensed fluid interface velocities

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-15262
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef FSI_MONOLITHICFLUIDSPLIT_H
#define FSI_MONOLITHICFLUIDSPLIT_H

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
  class MonolithicFluidSplit : public BlockMonolithic
  {

    friend class FSI::FSIResultTest;

  public:
    explicit MonolithicFluidSplit(const Epetra_Comm& comm,
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

    /// start a new time step
    virtual void PrepareTimeStep();

    /// recover Lagrange multiplier at the interface (i.e. condensed forces onto the fluid)
    /// needed for rhs in next time step
    virtual void RecoverLagrangeMultiplier();

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
    //! from the structure solution by a mapping across the interface.
    //! The condensed fluid degrees of freedom have to be recovered
    //! from the ale solution using a suitable displacement-velocity
    //! conversion.
    //!
    //! We are dealing with NOX here, so we get absolute values. x is the sum of
    //! all increments up to this point.
    //!
    //! \sa  ADAPTER::FluidFSI::DisplacementToVelocity()
    virtual void ExtractFieldVectors(
        Teuchos::RCP<const Epetra_Vector> x,   ///< composed vector that contains all field vectors
        Teuchos::RCP<const Epetra_Vector>& sx, ///< structural displacements
        Teuchos::RCP<const Epetra_Vector>& fx, ///< fluid velocities and pressure
        Teuchos::RCP<const Epetra_Vector>& ax  ///< ale displacements
    );


  private:

    /// build block vector from field vectors
    void SetupVector(
        Epetra_Vector &f,                      ///< composed vector that contains all field vectors
        Teuchos::RCP<const Epetra_Vector> sv,  ///< structural displacements
        Teuchos::RCP<const Epetra_Vector> fv,  ///< fluid velocities and pressure
        Teuchos::RCP<const Epetra_Vector> av,  ///< ale displacements
        const double fluidscale                ///< residual scaling for fluid
    );

    /// block system matrix
    Teuchos::RCP<OverlappingBlockMatrix> systemmatrix_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    Teuchos::RCP<UTILS::MatrixRowColTransform>  fggtransform_;
    Teuchos::RCP<UTILS::MatrixRowTransform>     fgitransform_;
    Teuchos::RCP<UTILS::MatrixColTransform>     figtransform_;
    Teuchos::RCP<UTILS::MatrixColTransform>     aigtransform_;
    Teuchos::RCP<UTILS::MatrixColTransform>     fmiitransform_;
    Teuchos::RCP<UTILS::MatrixRowColTransform>  fmgitransform_;

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

    /// @name Recovery of Lagrange multiplier at the end of each time step

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the fluid)
    //! evaluated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> lambda_;

    //! interface force \f$f_{\Gamma,i+1}^{F,n+1}\f$ onto the fluid at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const Epetra_Vector> fgcur_;

    //! interface force \f$f_{\Gamma,i}^{F,n+1}\f$ onto the fluid at previous NOX iteration \f$i\f$
    Teuchos::RCP<const Epetra_Vector> fgprev_;

    //! interface structure displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddginc_;

    //! inner fluid velocity increment \f$\Delta(\Delta u_{I,i+1}^{n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> duiinc_;

    //! interface displacement solution of the structure at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> disgprev_;

    //! inner velocity solution of fluid at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> soliprev_;

    //! interface velocity solution of the fluid at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> solgprev_;

    //! inner ALE displacement solution at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> solialeprev_;

    //! inner ALE displacement increment \f$\Delta(\Delta d_{I,i+1}^{G,n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddialeinc_;

    //! block \f$F_{\Gamma I,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fgicur_;

    //! block \f$F_{\Gamma I,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fgiprev_;

    //! block \f$F_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fggcur_;

    //! block \f$F_{\Gamma\Gamma,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fggprev_;

    //! block \f$F_{\Gamma I,i+1}^G\f$ of fluid shape derivatives matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmgicur_;

    //! block \f$F_{\Gamma I,i}^G\f$ of fluid shape derivatives matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmgiprev_;

    //! block \f$F_{\Gamma\Gamma,i+1}^G\f$ of fluid shape derivatives matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmggcur_;

    //! block \f$F_{\Gamma\Gamma,i}^G\f$ of fluid shape derivatives matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmggprev_;

    //@}
  };
}

#endif
