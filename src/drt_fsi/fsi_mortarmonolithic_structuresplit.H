/*----------------------------------------------------------------------------*/
/*!
\file fsi_mortarmonolithic_structuresplit.H

\maintainer Matthias Mayr

\brief Solve FSI problem with non-matching grids using a monolithic scheme
with condensed structure interface displacements

\level 1
*/
/*----------------------------------------------------------------------------*/

#ifndef FSI_MORTARMONOLITHICSTRUCTURESPLIT_H
#define FSI_MORTARMONOLITHICSTRUCTURESPLIT_H

#include "fsi_monolithic.H"
#include "../drt_inpar/inpar_fsi.H"

// forward declarations
namespace ADAPTER
{
  class Coupling;
  class CouplingMortar;
}

namespace  LINALG
{
  class BlockSparseMatrixBase;
}

namespace FSI
{
  class OverlappingBlockMatrix;

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

    In the sense of mortar coupling, structure split means that
    the structure field is chosen as slave field.
    Hence, the structural displacement interface degrees of freedom are
    condensed from the system along with the condensation of the Lagrange
    multiplier field, that is used to enforce the coupling conditions.

    The structural interface displacements are computed based on the
    fluid interface velocities. The conversion is done by
    ADAPTER::FluidFSI::VelocityToDisplacement().

    \sa MortarMonolithicStructureSplit
    \author tk
    \date 06/08
   */
  class MortarMonolithicStructureSplit : public BlockMonolithic
  {

    friend class FSI::FSIResultTest;

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

    //! Create #lambda_ and #lambdaold_
    virtual void SetLambda();

    //!Set #notsetup_ = true after redistribution
    virtual void SetNotSetup() { notsetup_ = true; return; }

    //! @name Apply current field state to system

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat);

    //@}

    /// the composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const;

    //! @name Methods for infnorm-scaling of the system

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(
        LINALG::BlockSparseMatrixBase& mat, ///< Jacobian matrix
        Epetra_Vector& b ///< right hand side
        );

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(
        LINALG::BlockSparseMatrixBase& mat, ///< Jacobian matrix
        Epetra_Vector& x, ///< solution vector
        Epetra_Vector& b ///< right hand side
        );

    //@}

    /// read restart
    virtual void ReadRestart(int step ///< step where we resatart from
        );

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
     *  = \left((a-b)\lambda^n + (b-a)\lambda^{n+1}\right)\left(d_\Gamma^{S,n+1}-d_\Gamma^{S,n}\right)
     *  \f]
     *  with the time interpolation factors a and b.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void CalculateInterfaceEnergyIncrement();

    /*! \brief Additional safety check of kinematic constraint during a single time step:
     *
     *  Constraint equation:
     *
     *  \f$D \mathbf{d}_{\Gamma}^{n+1} - D \mathbf{d}_{\Gamma}^{n} - \tau * M * \Delta \mathbf{u}_{\Gamma}^{n+1} - \Delta t M * \mathbf{u}_{\Gamma}^{n} \doteq \mathbf{0}\f$
     *
     *  with interface time integration factor
     *  \f$\tau = \begin{cases}\frac{\Delta t}{2} & \text {if }2^{nd}\text{ order}\\ \Delta t& \text {if }1^{st}\text{ order}\end{cases}\f$
     *
     *  Do this check only for safety reasons. Basically, the constraint is satisfied due to solving the
     *  condensed nonlinear system of equations. We expect really small violation norms.
     *
     *  \author mayr.mt \date 10/2012
     */
    virtual void CheckKinematicConstraint();

    /*! \brief Additional safety check of dynamic equilibrium during a single time step:
     *
     *  Dynamic equilibrium at the interface:
     *
     *  \f$M^{T} \mathbf{\lambda} - D^{T} \mathbf{\lambda} = \mathbf{0}\f$
     *
     *  Do this check only for safety reasons. Basically, the constraint is satisfied due to solving the
     *  condensed nonlinear system of equations. We expect really small violation norms.
     *
     *  \author mayr.mt \date 10/2012
     */
    virtual void CheckDynamicEquilibrium();

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
    virtual void CreateNodeOwnerRelationship(
        std::map<int,int>* nodeOwner,
        std::map<int,std::list<int> >* inverseNodeOwner,
        std::map<int, DRT::Node*>* structurenodesPtr,
        std::map<int, DRT::Node*>* fluidgnodesPtr,
        Teuchos::RCP<DRT::Discretization> structuredis,
        Teuchos::RCP<DRT::Discretization> fluiddis,
        const INPAR::FSI::Redistribute domain
        );

  protected:

    /// create the composed system matrix
    void CreateSystemMatrix();

    virtual void Update();

    virtual void Output();

    /// Write Lagrange multiplier
    virtual void OutputLambda();

    /// setup of NOX convergence tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo>
    CreateStatusTest(Teuchos::ParameterList& nlParams, ///< parameter list
        Teuchos::RCP<NOX::Epetra::Group> grp ///< the NOX group
        );

    /*! \brief Extract the three field vectors from a given composed vector
     *
     *  The condensed ale degrees of freedom have to be recovered
     *  from the fluid solution using a suitable velocity-displacement
     *  conversion. After that, the ale solution increment is projected onto the
     *  structure in mortar style where a possible structural predictor has to
     *  be considered.
     *
     *  We are dealing with NOX here, so we get absolute values. x is the sum of
     *  all increments up to this point.
     *
     *  \sa  ADAPTER::FluidFSI::VelocityToDisplacement()
     */
    virtual void ExtractFieldVectors(
        Teuchos::RCP<const Epetra_Vector> x, ///< composed vector that contains all field vectors
        Teuchos::RCP<const Epetra_Vector>& sx, ///< structural displacements
        Teuchos::RCP<const Epetra_Vector>& fx, ///< fluid velocities and pressure
        Teuchos::RCP<const Epetra_Vector>& ax ///< ale displacements
    );

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
    virtual void SetupDBCMapExtractor();

    /// setup RHS contributions based on single field residuals
    void SetupRHSResidual(Epetra_Vector& f);

    /// setup RHS contributions based on the Lagrange multiplier field
    void SetupRHSLambda(Epetra_Vector& f);

    /// setup RHS contributions based on terms for first nonlinear iteration
    void SetupRHSFirstiter(Epetra_Vector& f);

    /*! \brief Put three field vectors together to a single monolithic vector
     *
     *  If all input vectors contain interface DOFs, we first have to extract
     *  the inner portion of the slave vectors. Then, we just put them together.
     *  As usual, the ordering is: structure --  fluid -- ALE
     */
    void CombineFieldVectors(
          Epetra_Vector &v, ///< composed vector containing all field vectors
          Teuchos::RCP<const Epetra_Vector> sv, ///< structural DOFs
          Teuchos::RCP<const Epetra_Vector> fv, ///< fluid DOFs
          Teuchos::RCP<const Epetra_Vector> av, ///< ale DOFs
          bool fullvectors ///< all vectors contain all DOFs (true) or slave vectors only inner DOFs (false)
    );

    /// block system matrix
    Teuchos::RCP<OverlappingBlockMatrix> systemmatrix_;

    /// coupling of fluid and ale at the free surface
    Teuchos::RCP< ::ADAPTER::Coupling> fscoupfa_;

    /// coupling of structure and fluid at the interface
    Teuchos::RCP< ::ADAPTER::CouplingMortar> coupsfm_;

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

    /// @name Some quantities to recover the Lagrange multiplier at the end of each time step

    //! Lagrange multiplier \f$\lambda_\Gamma^n\f$ at the interface (ie condensed forces onto the structure)
    //! evaluated at old time step \f$t_n\f$ but needed for next time step \f$t_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> lambda_;

    //! Lagrange multiplier of previous time step
    Teuchos::RCP<Epetra_Vector> lambdaold_;

    //! inner structural displacement increment \f$\Delta(\Delta d_{I,i+1}^{n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> ddiinc_;

    //! interface fluid velocity increment \f$\Delta(\Delta u_{\Gamma,i+1}^{n+1})\f$ at current NOX iteration \f$i+1\f$
    Teuchos::RCP<Epetra_Vector> duginc_;

    //! inner displacement solution of the structure at previous NOX iteration
    Teuchos::RCP<const Epetra_Vector> disiprev_;

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

    //! summation of amount of artificial interface energy due to temporal discretization
    double energysum_;

    /// additional ale residual to avoid incremental ale errors
    Teuchos::RCP<Epetra_Vector> aleresidual_;

    /// preconditioned block Krylov or block Gauss-Seidel linear solver
    INPAR::FSI::LinearBlockSolver linearsolverstrategy_;

    /// ale movement relative to structure (none, slide_curr, slide_ref)
    INPAR::FSI::SlideALEProj aleproj_;

    bool notsetup_; ///< indicates if Setup has not been called yet

    Teuchos::RCP<FSI::UTILS::SlideAleUtils> slideale_; ///< Sliding Ale helper class

    Teuchos::RCP<Epetra_Vector> iprojdispinc_; ///<displacement of fluid side of the interface
    Teuchos::RCP<Epetra_Vector> iprojdisp_; ///<displacement of fluid side of the interface

  };
}

#endif
