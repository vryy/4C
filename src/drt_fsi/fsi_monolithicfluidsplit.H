/*----------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problem with matching grids using a monolithic scheme
with condensed fluid interface velocities

\maintainer Matthias Mayr

\level 1
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
}  // namespace ADAPTER

namespace NOX
{
  namespace FSI
  {
    class AdaptiveNewtonNormF;
    class Group;
  }  // namespace FSI
}  // namespace NOX

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
  }  // namespace UTILS
}  // namespace FSI

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
    explicit MonolithicFluidSplit(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /*! do the setup for the monolithic system


    1.) setup coupling; right now, we use matching meshes at the interface
    2.) create combined map
    3.) create block system matrix


    */
    virtual void SetupSystem();

    //! @name Apply current field state to system

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat);

    //@}

    /// the composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const;

    //! @name Methods for infnorm-scaling of the system

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b);

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(
        LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b);

    //@}

    /// start a new time step
    virtual void PrepareTimeStep();

    /*! \brief Recover Lagrange multiplier \f$\lambda_\Gamma\f$
     *
     *  Recover Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface at the
     *  end of each time step (i.e. condensed forces onto the structure) needed
     *  for rhs in next time step in order to guarantee temporal consistent
     *  exchange of coupling traction
     */
    virtual void RecoverLagrangeMultiplier();

    /*! \brief Compute spurious interface energy increment due to temporal discretization
     *
     *  Due to the temporal discretization, spurious energy \f$\Delta E_\Gamma^{n\rightarrow n+1}\f$
     *  might be produced at the interface. It can be computed as
     *  \f[
     *  \Delta E_\Gamma^{n\rightarrow n+1}
     *  = \left((a-b)\lambda^n +
     * (b-a)\lambda^{n+1}\right)\left(d_\Gamma^{S,n+1}-d_\Gamma^{S,n}\right) \f] with the time
     * interpolation factors a and b.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void CalculateInterfaceEnergyIncrement();

    /// Output routine accounting for Lagrange multiplier at the interface
    virtual void Output();

    /// Write Lagrange multiplier
    virtual void OutputLambda();

    //! take current results for converged and save for next time step
    virtual void Update();

    /// read restart data
    virtual void ReadRestart(int step);

    /// return Lagrange multiplier \f$\lambda_\Gamma\f$ at the interface
    Teuchos::RCP<Epetra_Vector> GetLambda() { return lambda_; };

    //! @name Time Adaptivity
    //@{

    /*! \brief Select \f$\Delta t_{min}\f$ of all proposed time step sizes based on error estimation
     *
     *  Depending on the chosen method (fluid or structure split), only 3 of the
     *  6 available norms are useful. Each of these three norms delivers a new
     *  time step size. Select the minimum of these three as the new time step size.
     *
     *  \author mayr.mt \date 08/2013
     */
    virtual double SelectDtErrorBased() const;

    /*! \brief Check whether time step is accepted or not
     *
     *  In case that the local truncation error is small enough, the time step is
     *  accepted.
     *
     *  \author mayr.mt \date 08/2013
     */
    virtual bool SetAccepted() const;

    //@}

    /*! \brief Find future / desired owner for each node at the interface
     *
     *  The relation is saved in the map \c nodeOwner as node -- owner.
     *
     *  In \c inverseNodeOwner the same information is contained in the form
     *  owner -- nodes.
     *
     *  The maps are built for interface nodes of the domain \c domain, where
     *  domain = {fluid, structure}.
     */
    virtual void CreateNodeOwnerRelationship(std::map<int, int>* nodeOwner,
        std::map<int, std::list<int>>* inverseNodeOwner, std::map<int, DRT::Node*>* fluidnodesPtr,
        std::map<int, DRT::Node*>* structuregnodesPtr,
        Teuchos::RCP<DRT::Discretization> structuredis, Teuchos::RCP<DRT::Discretization> fluiddis,
        const INPAR::FSI::Redistribute domain)
    {
      dserror("Not implemented!");
    }

   protected:
    /// create the composed system matrix
    void CreateSystemMatrix();

    /// setup of NOX convergence tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo> CreateStatusTest(
        Teuchos::ParameterList& nlParams, Teuchos::RCP<NOX::Epetra::Group> grp);

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
        Teuchos::RCP<const Epetra_Vector> x,    ///< composed vector that contains all field vectors
        Teuchos::RCP<const Epetra_Vector>& sx,  ///< structural displacements
        Teuchos::RCP<const Epetra_Vector>& fx,  ///< fluid velocities and pressure
        Teuchos::RCP<const Epetra_Vector>& ax   ///< ale displacements
    );

    /*! \brief Create the combined DOF row map for the FSI problem
     *
     *  Combine the DOF row maps of structure, fluid and ALE to an global FSI
     *  DOF row map.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void CreateCombinedDofRowMap();

    /*! set the Lagrange multiplier (e.g. after restart, to be called from subclass)
     *  \author kruse
     *  \date 05/2014
     */
    virtual void SetLambda(Teuchos::RCP<Epetra_Vector> lambdanew)
    {
#ifdef DEBUG
      if (lambdanew == Teuchos::null || !lambdanew->Map().PointSameAs(lambda_->Map()))
        dserror("Map failure! Attempting to assign invalid vector to lambda_.");
#endif
      lambda_ = lambdanew;
    }

   private:
    /*! \brief Setup the Dirichlet map extractor
     *
     *  Create a map extractor #dbcmaps_ for the Dirichlet degrees of freedom
     *  for the entire FSI problem. This is done just by combining the
     *  condition maps and other maps from structure, fluid and ALE to a FSI-global
     *  condition map and other map.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void SetupDBCMapExtractor();

    /// setup RHS contributions based on single field residuals
    void SetupRHSResidual(Epetra_Vector& f);

    /// setup RHS contributions based on the Lagrange multiplier field
    void SetupRHSLambda(Epetra_Vector& f);

    /// setup RHS contributions based on terms for first nonlinear iteration
    void SetupRHSFirstiter(Epetra_Vector& f);

    //! Put three field vectors together to a monolithic vector
    //!
    //! If all input vectors contain interface DOFs, we first have to extract
    //! the inner portion of the slave vectors. Then, we just put them together.
    //! As usual, the ordering is: structure --  fluid -- ALE
    void CombineFieldVectors(Epetra_Vector& v,  ///< composed vector containing all field vectors
        Teuchos::RCP<const Epetra_Vector> sv,   ///< structural DOFs
        Teuchos::RCP<const Epetra_Vector> fv,   ///< fluid DOFs
        Teuchos::RCP<const Epetra_Vector> av,   ///< ale DOFs
        bool fullvectors  ///< all vectors contain all DOFs (true) or slave vectors only inner DOFs
                          ///< (false)
    );

    /// block system matrix
    Teuchos::RCP<OverlappingBlockMatrix> systemmatrix_;

    /// @name Matrix block transform objects
    /// Handle row and column map exchange for matrix blocks

    Teuchos::RCP<UTILS::MatrixRowColTransform> fggtransform_;
    Teuchos::RCP<UTILS::MatrixRowTransform> fgitransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> figtransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> aigtransform_;
    Teuchos::RCP<UTILS::MatrixColTransform> fmiitransform_;
    Teuchos::RCP<UTILS::MatrixRowColTransform> fmgitransform_;
    Teuchos::RCP<UTILS::MatrixRowColTransform> fmggtransform_;

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

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the
    //! fluid) evaluated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> lambda_;

    //! Lagrange multiplier of previous time step
    Teuchos::RCP<Epetra_Vector> lambdaold_;

    //! interface structure displacement increment \f$\Delta(\Delta d_{\Gamma,i+1}^{n+1})\f$ at
    //! current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddginc_;

    //! inner fluid velocity increment \f$\Delta(\Delta u_{I,i+1}^{n+1})\f$ at current NOX iteration
    //! \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> duiinc_;

    //! interface displacement solution of the structure at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> disgprev_;

    //! inner velocity solution of fluid at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> soliprev_;

    //! interface velocity solution of the fluid at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> solgprev_;

    //! inner ALE displacement solution at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> solialeprev_;

    //! inner ALE displacement increment \f$\Delta(\Delta d_{I,i+1}^{G,n+1})\f$ at current NOX
    //! iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddialeinc_;

    //! block \f$F_{\Gamma I,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fgicur_;

    //! block \f$F_{\Gamma I,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fgiprev_;

    //! block \f$F_{\Gamma\Gamma,i+1}\f$ of fluid matrix at current NOX iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fggcur_;

    //! block \f$F_{\Gamma\Gamma,i}\f$ of fluid matrix at previous NOX iteration \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fggprev_;

    //! block \f$F_{\Gamma I,i+1}^G\f$ of fluid shape derivatives matrix at current NOX iteration
    //! \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmgicur_;

    //! block \f$F_{\Gamma I,i}^G\f$ of fluid shape derivatives matrix at previous NOX iteration
    //! \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmgiprev_;

    //! block \f$F_{\Gamma\Gamma,i+1}^G\f$ of fluid shape derivatives matrix at current NOX
    //! iteration \f$i+1\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmggcur_;

    //! block \f$F_{\Gamma\Gamma,i}^G\f$ of fluid shape derivatives matrix at previous NOX iteration
    //! \f$i\f$
    Teuchos::RCP<const LINALG::SparseMatrix> fmggprev_;

    //@}

    //! summation of amount of artificial interface energy due to temporal discretization
    double energysum_;
  };
}  // namespace FSI

#endif
