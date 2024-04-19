/*----------------------------------------------------------------------*/
/*! \file
\brief Volume-coupled FSI (base class)

\level 3

*/
/*----------------------------------------------------------------------*/



#ifndef FOUR_C_FSI_LUNGMONOLITHIC_HPP
#define FOUR_C_FSI_LUNGMONOLITHIC_HPP

#include "baci_config.hpp"

#include "baci_coupling_adapter.hpp"
#include "baci_fsi_monolithic.hpp"
#include "baci_inpar_fsi.hpp"

#include <fstream>

FOUR_C_NAMESPACE_OPEN

// forward declarations

namespace CONSTRAINTS
{
  class ConstraintDofSet;
}

namespace CORE::LINALG
{
  class BlockSparseMatrixBase;
}

namespace ADAPTER
{
  class Coupling;
}

namespace FSI
{
  // forward declaration
  class OverlappingBlockMatrix;

  /// monolithic FSI algorithm with overlapping interface equations
  /// for simulation of a specific class of bio problems (FSI airway
  /// model with attached balloon built of lung parenchyma)
  class LungMonolithic : public BlockMonolithic
  {
   public:
    explicit LungMonolithic(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    /*! do the setup for the monolithic system
      1.) setup coupling; right now, we use matching meshes at the interface
      2.) create combined map
      3.) create block system matrix
    */
    void SetupSystem() override = 0;

    /// some general setup stuff necessary for both fluid and
    /// structure split
    void GeneralSetup();

    /// setup composed system matrix from field solvers
    void SetupSystemMatrix(CORE::LINALG::BlockSparseMatrixBase& mat) override = 0;

    //@}

    /// Evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    void Evaluate(
        Teuchos::RCP<const Epetra_Vector> step_increment  ///< increment between time step n and n+1
        ) override;

    /// Extract initial guess from fields
    void InitialGuess(Teuchos::RCP<Epetra_Vector> ig) override = 0;

    /// the composed system matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> SystemMatrix() const override;

    //! @name Methods for infnorm-scaling of the system

    /// apply infnorm scaling to linear block system
    void ScaleSystem(CORE::LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b) override;

    /// undo infnorm scaling from scaled solution
    void UnscaleSolution(
        CORE::LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b) override;

    //@}

    /// Update everything
    void Update() override;

    CORE::ADAPTER::Coupling& StructureAleOutflowCoupling() { return *coupsaout_; }

    void ReadRestart(int step) override;

    //! @name Time Adaptivity
    //@{

    /*! \brief Select \f$\Delta t_{min}\f$ of all proposed time step sizes based on error estimation
     *
     *  Depending on the chosen method (fluid or structure split), only 3 of the
     *  6 available norms are useful. Each of these three norms delivers a new
     *  time step size. Select the minimum of these three as the new time step size.
     */
    double SelectDtErrorBased() const override
    {
      FOUR_C_THROW("SelectDtErrorBased() not implemented, yet!");
      return 0.0;
    }

    /*! \brief Check whether time step is accepted or not
     *
     *  In case that the local truncation error is small enough, the time step is
     *  accepted.
     */
    bool SetAccepted() const override
    {
      FOUR_C_THROW("SetAccepted() not implemented, yet!");
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
    void CreateNodeOwnerRelationship(std::map<int, int>* nodeOwner,
        std::map<int, std::list<int>>* inverseNodeOwner, std::map<int, DRT::Node*>* fluidnodesPtr,
        std::map<int, DRT::Node*>* structuregnodesPtr,
        Teuchos::RCP<DRT::Discretization> structuredis, Teuchos::RCP<DRT::Discretization> fluiddis,
        const INPAR::FSI::Redistribute domain) override
    {
      FOUR_C_THROW("Not implemented, yet.");
    }

   protected:
    /// create the composed system matrix
    void CreateSystemMatrix(bool structuresplit);

    void Output() override;

    /// start a new time step
    void PrepareTimeStep() override;

    /// transfer helper
    Teuchos::RCP<Epetra_Vector> StructToAleOutflow(Teuchos::RCP<Epetra_Vector> iv) const;

    /// setup solver for global block system
    Teuchos::RCP<::NOX::Epetra::LinearSystem> CreateLinearSystem(Teuchos::ParameterList& nlParams,
        ::NOX::Epetra::Vector& noxSoln, Teuchos::RCP<::NOX::Utils> utils) override;

    /// setup of NOX convergence tests
    Teuchos::RCP<::NOX::StatusTest::Combo> CreateStatusTest(
        Teuchos::ParameterList& nlParams, Teuchos::RCP<::NOX::Epetra::Group> grp) override;

    /// extract the three field vectors from a given composed vector
    /*!
      We are dealing with NOX here, so we get absolute values. x is the sum of
      all increments up to this point.

      \param x  (i) composed vector that contains all field vectors
      \param sx (o) structural displacements
      \param fx (o) fluid velocities and pressure
      \param ax (o) ale displacements
    */
    void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
        Teuchos::RCP<const Epetra_Vector>& ax) override = 0;

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
    Teuchos::RCP<CORE::ADAPTER::Coupling> icoupfa_;

    /// additional coupling of structure and ale fields at airway outflow
    Teuchos::RCP<CORE::ADAPTER::Coupling> coupsaout_;

    /// additional coupling of structure and ale/fluid fields at airway outflow
    Teuchos::RCP<CORE::ADAPTER::Coupling> coupfsout_;

    /// fluid and ale coupling at airway outflow
    Teuchos::RCP<CORE::ADAPTER::Coupling> coupfaout_;

    ///@}

    /// @name infnorm scaling

    Teuchos::RCP<Epetra_Vector> srowsum_;
    Teuchos::RCP<Epetra_Vector> scolsum_;
    Teuchos::RCP<Epetra_Vector> arowsum_;
    Teuchos::RCP<Epetra_Vector> acolsum_;

    //@}

    ///@}

    /// @lung fluid-structure volume constraints

    Teuchos::RCP<Epetra_Vector> LagrMultVec_;     ///< lagrange multipliers
    Teuchos::RCP<Epetra_Vector> LagrMultVecOld_;  ///< lagrange multipliers of last time step
    Teuchos::RCP<Epetra_Vector>
        IncLagrMultVec_;  ///< nonlinear iteration increment vector of lagrange multipliers

    Teuchos::RCP<CONSTRAINTS::ConstraintDofSet>
        ConstrDofSet_;                       ///< degrees of freedom of lagrange multipliers
    int OffsetID_;                           ///< smallest constraint boundary condition ID
    int NumConstrID_;                        ///< number of volume constraint boundary conditions
    int NumPresConstrID_;                    ///< number of pressure constraint boundary conditions
    Teuchos::RCP<Epetra_Map> ConstrMap_;     ///< unique map of constraint values
    Teuchos::RCP<Epetra_Map> RedConstrMap_;  ///< fully redundant map of constraint values
    Teuchos::RCP<Epetra_Export>
        ConstrImport_;  ///< importer fully redundant <-> unique map of constraint values

    Teuchos::RCP<Epetra_Vector> SignVolsRed_;  ///< signs of volumes
    Teuchos::RCP<Epetra_Vector>
        OldVols_;  ///< volumes of last time step needed for rhs of constraint equations
    Teuchos::RCP<Epetra_Vector>
        CurrVols_;  ///< current volumes needed for rhs of constraint equations
    Teuchos::RCP<Epetra_Vector>
        OldFlowRates_;  ///< flow rates of last time step needed for rhs of constraint equations
    Teuchos::RCP<Epetra_Vector>
        CurrFlowRates_;  ///< current flow rates needed for rhs of constraint equations

    Teuchos::RCP<Epetra_Vector>
        dVfluid_;  ///< current change in fluid volumes (for output purposes only)
    Teuchos::RCP<Epetra_Vector>
        dVstruct_;  ///< current change in structure volumes (for output purposes only)

    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>
        AddStructConstrMatrix_;  ///< matrix containing all structure constraint related stuff

    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>
        AddFluidShapeDerivMatrix_;  ///< additional constraint portion on block (1,2)
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        FluidConstrMatrix_;  ///< rectangular fluid matrix associated with constraints K_fl
    Teuchos::RCP<CORE::LINALG::SparseMatrix>
        ConstrFluidMatrix_;  ///< rectangular fluid matrix associated with constraints K_lf

    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>
        AleConstrMatrix_;  ///< rectangular ale matrix associated with constraints K_al
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>
        ConstrAleMatrix_;  ///< rectangular ale matrix associated with constraints K_la

    Teuchos::RCP<Epetra_Vector> AddStructRHS_;  ///< additional constraint portion on structural rhs
    Teuchos::RCP<Epetra_Vector> AddFluidRHS_;   ///< additional constraint portion on fluid rhs
    Teuchos::RCP<Epetra_Vector> ConstrRHS_;     ///< rhs of constraint equations

    double theta_;  ///< parameter for integrating fluid volumes (one-step theta scheme)

    //@}

    /// preconditioned block Krylov or block Gauss-Seidel linear solver
    INPAR::FSI::LinearBlockSolver linearsolverstrategy_;

    /// output of changes in volumes in text file
    std::ofstream outfluiddvol_;
    std::ofstream outstructdvol_;
    std::ofstream outstructabsvol_;

   private:
    /*! \brief Create the combined DOF row map for the FSI problem
     *
     *  Combine the DOF row maps of structure, fluid and ALE to an global FSI
     *  DOF row map.
     */
    void CreateCombinedDofRowMap() override = 0;

    /*! \brief Setup the Dirichlet map extractor
     *
     *  Create a map extractor #dbcmaps_ for the Dirichlet degrees of freedom
     *  for the entire FSI problem. This is done just by combining the
     *  condition maps and other maps from structure, fluid and ALE to a FSI-global
     *  condition map and other map.
     */
    void SetupDBCMapExtractor() override = 0;

    /// setup RHS contributions based on single field residuals
    void SetupRHSResidual(Epetra_Vector& f) override = 0;

    /// setup RHS contributions based on the Lagrange multiplier field
    void SetupRHSLambda(Epetra_Vector& f) override = 0;

    /// setup RHS contributions based on terms for first nonlinear iteration
    void SetupRHSFirstiter(Epetra_Vector& f) override = 0;
  };
}  // namespace FSI

FOUR_C_NAMESPACE_CLOSE

#endif
