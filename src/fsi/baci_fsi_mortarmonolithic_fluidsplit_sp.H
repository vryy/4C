/*----------------------------------------------------------------------------*/
/*! \file

\brief Solve FSI problem with non-matching grids using a monolithic scheme
in saddle-point formulation with Lagrange multipliers discretized on the fluid interface

\level 2
*/

/*----------------------------------------------------------------------------*/

#ifndef BACI_FSI_MORTARMONOLITHIC_FLUIDSPLIT_SP_H
#define BACI_FSI_MORTARMONOLITHIC_FLUIDSPLIT_SP_H

#include "baci_fsi_monolithic.H"
#include "baci_inpar_fsi.H"

// forward declarations
namespace CORE::ADAPTER
{
  class Coupling;
  class CouplingMortar;
}  // namespace CORE::ADAPTER

class Epetra_Comm;

namespace CORE::LINALG
{
  class BlockSparseMatrixBase;
  class MatrixColTransform;
}  // namespace CORE::LINALG

namespace FSI
{
  class OverlappingBlockMatrix;

  namespace UTILS
  {
    class SlideAleUtils;
  }  // namespace UTILS
}  // namespace FSI

namespace NOX
{
  namespace Epetra
  {
    class Group;
  }

  namespace StatusTest
  {
    class Combo;
  }
}  // namespace NOX

namespace Teuchos
{
  class ParameterList;
}

namespace FSI
{
  class MortarMonolithicFluidSplitSaddlePoint : public BlockMonolithic
  {
    friend class FSI::FSIResultTest;

   public:
    explicit MortarMonolithicFluidSplitSaddlePoint(
        const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams);

    void SetupSystem() final;

    //! @name Apply current field state to system

    /// setup composed system matrix from field solvers
    void SetupSystemMatrix(CORE::LINALG::BlockSparseMatrixBase& mat) final;

    //@}

    /// the composed system matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> SystemMatrix() const override;

    /// read restart
    void ReadRestart(int step) final;

    //! @name Time Adaptivity
    //@{

    /*! \brief Select \f$\Delta t_{min}\f$ of all proposed time step sizes
     *         based on error estimation
     *
     *  Depending on the chosen method (fluid or structure split), only 3 of the
     *  6 available norms are useful. Each of these three norms delivers a new
     *  time step size. Select the minimum of these three as the new time step
     *  size.
     */
    double SelectDtErrorBased() const final;

    /*! \brief Check whether time step is accepted or not
     *
     *  In case that the local truncation error is small enough, the time step
     *  is accepted.
     */
    bool SetAccepted() const final;

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
        const INPAR::FSI::Redistribute domain) final;

   protected:
    void CreateSystemMatrix();

    Teuchos::RCP<NOX::StatusTest::Combo> CreateStatusTest(
        Teuchos::ParameterList& nlParams, Teuchos::RCP<NOX::Epetra::Group> grp) final;

    void Update() final;

    void Output() final;

    /// Write Lagrange multiplier
    void OutputLambda() final;

    /*!
    @copydoc FSI::Monolithic::ExtractFieldVectors

    Since this is a saddle-point formulation, we also need to extract the Lagrange multipliers
    #lag_mult_ from the monolithic solution vector.
    */
    void ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx,
        Teuchos::RCP<const Epetra_Vector>& ax, Teuchos::RCP<const Epetra_Vector>& lagx);


    /*!
    @copydoc FSI::Monolithic::Evaluate

    Since this is a saddle-point formulation, we also need to evaluate the Lagrange multipliers
    #lag_mult_.
    */
    void Evaluate(
        Teuchos::RCP<const Epetra_Vector> stepinc  ///< increment between time step n and n+1
        ) override;

    void InitialGuess(Teuchos::RCP<Epetra_Vector> initial_guess) final;

   private:
    /*! \brief Create the combined DOF row map for the FSI problem
     *
     *  Combine the DOF row maps of structure, fluid, ALE and Lagrange multipliers to an global FSI
     *  DOF row map.
     */
    void CreateCombinedDofRowMap() final;

    /*! \brief Create the DOF row map for lagrange
     *
     *  Create the DOF row map for lagrange multiplier based on the fluid interface field and
     *  last GID of the ALE field DOF row map
     */
    virtual void CreateLagrangeMultiplierDofRowMap();

    virtual void CombineFieldVectors(
        Epetra_Vector& f,  ///< composed vector containing all field vectors
        Teuchos::RCP<const Epetra_Vector> solid_vector,     ///< structural DOFs
        Teuchos::RCP<const Epetra_Vector> fluid_vector,     ///< fluid DOFs
        Teuchos::RCP<const Epetra_Vector> ale_vector,       ///< ale DOFs
        Teuchos::RCP<const Epetra_Vector> lag_mult_vector,  /// < lagrange multiplier
        bool fullvectors);

    /*! \brief Setup the Dirichlet map extractor
     *
     *  Create a map extractor #dbcmaps_ for the Dirichlet degrees of freedom
     *  for the entire FSI problem. This is done just by combining the
     *  condition maps and other maps from structure, fluid and ALE to a
     *  FSI-global condition map and other map.
     */
    void SetupDBCMapExtractor() final;

    /// setup RHS contributions based on single field residuals
    void SetupRHSResidual(Epetra_Vector& f) final;

    /// setup RHS contributions based on the Lagrange multiplier field
    void SetupRHSLambda(Epetra_Vector& f) final;

    /// setup RHS contributions based on terms for first nonlinear iteration
    void SetupRHSFirstiter(Epetra_Vector& f) final;

    //! Create #lag_mult_
    virtual void SetLagMult();

    //! Set #notsetup_ = true after redistribution
    void SetNotSetup() override { notsetup_ = true; }

    /*! block system matrix
     *  System matrix has a 6x6-block structure corresponding to the vector of unknowns
     *
     *  \f$\Delta x^T = [\Delta d_I^{S,n+1}~\Delta d_\Gamma^{S,n+1}~\Delta u_I^{F,n+1}~\Delta
     * u_\Gamma^{F,n+1}~\Delta d_I^{G,n+1}~\Lambda^{n+1}]\f$.
     *
     * As for the code, we have a 4x4 system.
     */
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> systemmatrix_;

    /// communicator
    const Epetra_Comm& comm_;

    /// @name Matrix block transform objects to handle row and column map exchange for matrix blocks

    /// Coupling of structure and fluid at the interface
    Teuchos::RCP<CORE::ADAPTER::CouplingMortar> coupling_solid_fluid_mortar_;

    /// Helper variable for the transformation of aleunknowns onto the slave side
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> ale_inner_interf_transform_;

    /// Helper variable for the transformation of fluid unknowns onto the slave side
    Teuchos::RCP<CORE::LINALG::MatrixColTransform> fluid_mesh_inner_inner_transform_;

    ///@}

    /// additional ale residual to avoid incremental ale errors
    Teuchos::RCP<Epetra_Vector> aleresidual_;

    //! DOF map of Lagrange multiplier unknowns
    Teuchos::RCP<const Epetra_Map> lag_mult_dof_map_;

    //! Lagrange multiplier
    Teuchos::RCP<Epetra_Vector> lag_mult_;

    //! Lagrange multiplier from previous time step
    Teuchos::RCP<Epetra_Vector> lag_mult_old_;

    //! Flag to indicate if Setup has not been called yet
    bool notsetup_;

  };  // class MortarMonolithicFluidSplitSaddlePoint
}  // namespace FSI

#endif  // FSI_MORTARMONOLITHIC_FLUIDSPLIT_SP_H