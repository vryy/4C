/*----------------------------------------------------------------------*/
/*!

\brief Solve FSI problem with constraints

\level 2

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------*/

#ifndef FSI_CONSTR_MONOLITHIC_H
#define FSI_CONSTR_MONOLITHIC_H

#include "fsi_constr_overlapprec.H"
#include "fsi_monolithic.H"

// forward declarations
namespace UTILS
{
  class ConstrManager;
}

namespace ADAPTER
{
  class Coupling;
}

namespace FSI
{
  /// monolithic FSI algorithm with overlapping interface equations
  /// for simulation of a algebraically constrained structure field
  /*!
    \author tk
    \date 02/10
   */
  class ConstrMonolithic : public BlockMonolithic
  {
   public:
    explicit ConstrMonolithic(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /*! do the setup for the monolithic system
        1.) setup coupling; right now, we use matching meshes at the interface
    2.) create combined map
    3.) create block system matrix
        */
    virtual void SetupSystem() = 0;

    /// some general setup stuff necessary for both fluid and
    /// structure split
    void GeneralSetup();

    /// setup composed system matrix from field solvers
    virtual void SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat) = 0;

    //@}

    /// Evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc  ///< increment between time step n and n+1
    );

    /// Extract initial guess from fields
    virtual void InitialGuess(Teuchos::RCP<Epetra_Vector> ig) = 0;

    /// the composed system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> SystemMatrix() const { return systemmatrix_; }

    //! @name Methods for infnorm-scaling of the system

    /// apply infnorm scaling to linear block system
    virtual void ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b);

    /// undo infnorm scaling from scaled solution
    virtual void UnscaleSolution(
        LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b);

    //@}

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
    virtual double SelectDtErrorBased() const
    {
      dserror("SelectDtErrorBased() not implemented, yet!");
      return 0.0;
    }

    /*! \brief Check whether time step is accepted or not
     *
     *  In case that the local truncation error is small enough, the time step is
     *  accepted.
     *
     *  \author mayr.mt \date 08/2013
     */
    virtual bool SetAccepted() const
    {
      dserror("SetAccepted() not implemented, yet!");
      return false;
    }

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
      dserror("Not implemented, yet.");
    }

   protected:
    /// create the composed system matrix
    void CreateSystemMatrix(bool structuresplit);

    /// setup solver for global block system
    Teuchos::RCP<NOX::Epetra::LinearSystem> CreateLinearSystem(Teuchos::ParameterList& nlParams,
        NOX::Epetra::Vector& noxSoln, Teuchos::RCP<NOX::Utils> utils);

    /// setup of NOX convergence tests
    virtual Teuchos::RCP<NOX::StatusTest::Combo> CreateStatusTest(
        Teuchos::ParameterList& nlParams, Teuchos::RCP<NOX::Epetra::Group> grp);

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
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
        Teuchos::RCP<const Epetra_Vector>& ax) = 0;

    /// build block vector from field vectors
    virtual void SetupVector(Epetra_Vector& f,
        Teuchos::RCP<const Epetra_Vector> sv,  ///< structure vector
        Teuchos::RCP<const Epetra_Vector> fv,  ///< fluid vector
        Teuchos::RCP<const Epetra_Vector> av,  ///< ale vector
        Teuchos::RCP<const Epetra_Vector> cv,  ///< constraint vector
        double fluidscale) = 0;                ///< scaling

    /// block system matrix
    Teuchos::RCP<OverlappingBlockMatrix> systemmatrix_;

    /// ALE residual
    Teuchos::RCP<Epetra_Vector> aleresidual_;

    /// restart information
    int writerestartevery_;

    /// coupling of fluid and ale (interface only)
    Teuchos::RCP<ADAPTER::Coupling> icoupfa_;

    /// additional coupling of structure and ale fields at airway outflow
    Teuchos::RCP<ADAPTER::Coupling> coupsaout_;

    /// additional coupling of structure and ale/fluid fields at airway outflow
    Teuchos::RCP<ADAPTER::Coupling> coupfsout_;

    /// fluid and ale coupling at airway outflow
    Teuchos::RCP<ADAPTER::Coupling> coupfaout_;

    ///@}

    /// @name infnorm scaling

    Teuchos::RCP<Epetra_Vector> srowsum_;
    Teuchos::RCP<Epetra_Vector> scolsum_;
    Teuchos::RCP<Epetra_Vector> arowsum_;
    Teuchos::RCP<Epetra_Vector> acolsum_;

    //@}

    ///@}

    /// @ information about constraints are taken from structure ConstraintManager

    const Teuchos::RCP<::UTILS::ConstrManager> conman_;  ///< constraint manager

    //@}

    /// preconditioned block Krylov or block Gauss-Seidel linear solver
    INPAR::FSI::LinearBlockSolver linearsolverstrategy_;


   private:
    /*! \brief Create the combined DOF row map for the FSI problem
     *
     *  Combine the DOF row maps of structure, fluid and ALE to an global FSI
     *  DOF row map.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void CreateCombinedDofRowMap() = 0;

    /*! \brief Setup the Dirichlet map extractor
     *
     *  Create a map extractor #dbcmaps_ for the Dirichlet degrees of freedom
     *  for the entire FSI problem. This is done just by combining the
     *  condition maps and other maps from structure, fluid and ALE to a FSI-global
     *  condition map and other map.
     *
     *  \author mayr.mt \date 05/2014
     */
    virtual void SetupDBCMapExtractor() = 0;

    /// setup RHS contributions based on single field residuals
    virtual void SetupRHSResidual(Epetra_Vector& f) = 0;

    /// setup RHS contributions based on the Lagrange multiplier field
    virtual void SetupRHSLambda(Epetra_Vector& f) = 0;

    /// setup RHS contributions based on terms for first nonlinear iteration
    virtual void SetupRHSFirstiter(Epetra_Vector& f) = 0;
  };
}  // namespace FSI

#endif
