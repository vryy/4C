/*----------------------------------------------------------------------*/
/*!
\file xfluidfluid.H

\brief Control routine for fluid-fluid (in)stationary solvers with XFEM

\level 2

<pre>
\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
 */
/*----------------------------------------------------------------------*/

#ifndef XFLUIDFLUID_H_
#define XFLUIDFLUID_H_


#include "xfluid.H"
#include "xfluidfluid_state.H"


namespace XFEM
{
  class MeshProjector;
  class MeshCouplingFluidFluid;
}  // namespace XFEM

namespace FLD
{
  namespace UTILS
  {
    class XFluidFluidMapExtractor;
  }

  /*!
    Class can handle a fluid described on a XFEM discretization and an embedded
    standard FEM discretization

    \author kruse
    \date 01/15
  */
  class XFluidFluid : public XFluid
  {
    friend class XFluidResultTest;

   public:
    //! ctor
    XFluidFluid(const Teuchos::RCP<FLD::FluidImplicitTimeInt>& embedded_fluid,  ///< embedded fluid
        const Teuchos::RCP<DRT::Discretization>& xfluiddis,  ///< background fluid discretization
        const Teuchos::RCP<LINALG::Solver>& solver,          ///< fluid solver
        const Teuchos::RCP<Teuchos::ParameterList>& params,  ///< xfluid params
        bool ale_xfluid = false,  ///< background (XFEM) fluid in ALE-formulation
        bool ale_fluid = false    ///< embedded fluid in ALE-formulation
    );

    //! ctor for multiple mesh coupling
    XFluidFluid(const Teuchos::RCP<FLD::FluidImplicitTimeInt>& embedded_fluid,  ///< embedded fluid
        const Teuchos::RCP<DRT::Discretization>& xfluiddis,  ///< background fluid discretization
        const Teuchos::RCP<DRT::Discretization>&
            soliddis,                                ///< structure discretization to couple with
        const Teuchos::RCP<LINALG::Solver>& solver,  ///< fluid solver
        const Teuchos::RCP<Teuchos::ParameterList>& params,  ///< xfluid params
        bool ale_xfluid = false,  ///< background (XFEM) fluid in ALE-formulation
        bool ale_fluid = false    ///< embedded fluid in ALE-formulation
    );

    virtual ~XFluidFluid();

    /// initialization
    virtual void Init() { Init(true); }
    virtual void Init(bool createinitialstate);

    virtual void CreateInitialState();

    /// don't keep a sparse matrix underneath
    void UseBlockMatrix(bool splitmatrix = true);

    /// set fluid-fluid coupling specific parameters
    void SetXFluidFluidParams();

    /// set initial flow field for fluid domains
    void SetInitialFlowField(const INPAR::FLUID::InitialField initfield, const int startfuncno);

    /// set fluid-fluid interface fixed for current time step
    void SetInterfaceFixed();

    /// free fluid-fluid interface
    void SetInterfaceFree();

    /// setup the variables to do a new time step
    void PrepareTimeStep();

    /// get initial guess for both fluids
    Teuchos::RCP<const Epetra_Vector> InitialGuess();

    /// prepare solution (cut happens here)
    void PrepareXFEMSolve();

    /// Monolithic FSI needs to access the linear fluid problem.
    void Evaluate(Teuchos::RCP<const Epetra_Vector>
            stepinc  ///< solution increment between time step n and n+1
    );

    /// Update the solution after convergence of the nonlinear
    /// iteration. Current solution becomes old solution of next timestep.
    void TimeUpdate();

    /*!
     * \brief set underlying dof-maps for new shape derivatives matrix
     * (required in the context of monolithic XFFSI)
     * \param (in) fsiextractor map extractor to distinguish between FSI/non-FSI dof
     * in the merged fluid-fluid-dof-map
     * \param (in) condelements map of conditioned elements
     */
    void PrepareShapeDerivatives(const Teuchos::RCP<const LINALG::MultiMapExtractor> fsiextractor,
        const Teuchos::RCP<std::set<int>> condelements);

    /*!
     * \brief return fluid-fluid system matrix as block matrix
     * (required in the context of monolithic fluidsplit-XFFSI)
     * \param (in) innermap map of inner embedded and background fluid dof
     * \param (in) condmap map of fsi interface dof
     * \return coupled fluid-fluid block system matrix
     */
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix(
        Teuchos::RCP<Epetra_Map> innermap, Teuchos::RCP<Epetra_Map> condmap);

    Teuchos::RCP<const Epetra_Vector> RHS() { return xff_state_->xffluidresidual_; }
    Teuchos::RCP<const Epetra_Vector> Velnp() { return xff_state_->xffluidvelnp_; }
    Teuchos::RCP<Epetra_Vector> WriteAccessVelnp() { return xff_state_->xffluidvelnp_; }
    Teuchos::RCP<const Epetra_Vector> Veln() { return xff_state_->xffluidveln_; }
    Teuchos::RCP<const Epetra_Vector> Stepinc() { return stepinc_; }

    Teuchos::RCP<Epetra_Vector> WriteAccessDispOldState() { return dispnpoldstate_; }

    /// get merged fluid dofmap
    Teuchos::RCP<const Epetra_Map> DofRowMap()
    {
      //    return Teuchos::rcp((xff_state_->xffluiddofrowmap_).get(), false); // TODO: does not
      //    work at the moment, no nox structure split gives segfault!
      return xff_state_->xffluiddofrowmap_;
    }

    /// get merged vel-pres-splitter
    Teuchos::RCP<LINALG::MapExtractor> VelPresSplitter()
    {
      return xff_state_->xffluidvelpressplitter_;
    }

    // get merged pressure dof-map
    Teuchos::RCP<const Epetra_Map> PressureRowMap();
    // get merged velocity dof-map
    Teuchos::RCP<const Epetra_Map> VelocityRowMap();

    Teuchos::RCP<LINALG::BlockSparseMatrixBase> ExtendedShapeDerivatives()
    {
      return extended_shapederivatives_;
    }

    /// get the combined fluid-fluid system matrix
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix()
    {
      return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(xff_state_->xffluidsysmat_);
    }

    /// get the combined fluid-fluid system matrix in block form (embedded and background fluid
    /// blocks)
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix()
    {
      return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(xff_state_->xffluidsysmat_);
    }

    /// get map extractor for background/embedded fluid
    Teuchos::RCP<FLD::UTILS::XFluidFluidMapExtractor> const& XFluidFluidMapExtractor()
    {
      return xff_state_->xffluidsplitter_;
    }

    //! get combined background and embedded fluid dirichlet map extractor
    Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
    {
      return xff_state_->xffluiddbcmaps_;
    }

    //! access to new merged fluid state container
    Teuchos::RCP<FLD::XFluidState> GetNewState();

    /// solve fluid field after ALE-mesh relaxation
    void UpdateMonolithicFluidSolution(const Teuchos::RCP<const Epetra_Map>& fsidofmap);

    /// interpolate embedded fluid state vectors from previous mesh displacement state
    void InterpolateEmbeddedStateVectors();

    /// create a result test
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    /// write output for both fluid discretizations
    void Output();

    /// evaluate errors compared to implemented analytical solutions
    Teuchos::RCP<std::vector<double>> EvaluateErrorComparedToAnalyticalSol();

   private:
    /// projection of history from other discretization - returns true if projection was successful
    /// for all nodes
    bool XTimint_ProjectFromEmbeddedDiscretization(
        const Teuchos::RCP<XFEM::XFluidTimeInt>& xfluid_timeint,  ///< xfluid time integration class
        std::vector<Teuchos::RCP<Epetra_Vector>>&
            newRowStateVectors,  ///< vectors to be reconstructed
        Teuchos::RCP<const Epetra_Vector>
            target_dispnp,     ///< displacement col - vector timestep n+1
        const bool screen_out  ///< screen output?
    );

    //! transfer dofs between Newton steps
    bool XTimint_DoIncrementStepTransfer(const bool screen_out, const bool firstcall_in_timestep);

    //! create a new merged state container
    void CreateState();

    /// call the loop over elements to assemble volume and interface integrals
    void AssembleMatAndRHS(int itnum  ///< iteration number
    );

    /// Update velocity and pressure by increment
    void UpdateByIncrement();

    /// add additional EOS terms to interface-contributing embedded fluid elements
    void AddEosPresStabToEmbLayer();

    //! @name XFF time integration

    /// embedded fluid field
    Teuchos::RCP<FLD::FluidImplicitTimeInt> embedded_fluid_;

    /// pointer to fluid-fluid mesh coupling object
    Teuchos::RCP<XFEM::MeshCouplingFluidFluid> mc_xff_;

    /// for XFEM time integration by projection from embedded mesh
    Teuchos::RCP<XFEM::MeshProjector> projector_;

    /// whether the embedded fluid is an ALE-fluid
    bool ale_embfluid_;

    /// merged fluid state (cut-dependent) at time n+1
    Teuchos::RCP<FLD::XFluidFluidState> xff_state_;

    /// vector with Newton increments (used for monolithic fully newton fsi approach)
    Teuchos::RCP<Epetra_Vector> stepinc_;

    /// ALE-displacements of previous step
    Teuchos::RCP<Epetra_Vector> dispnpoldstate_;

    /// shape derivatives matrix (linearization with respect to mesh motion),
    /// including background fluid dof (set to zero)
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> extended_shapederivatives_;

    /// flag, that indicates active shape derivatives
    bool active_shapederivatives_;

    /// flag, that indicates evaluation of pressure EOS-terms on outer embedded elements
    bool xff_eos_pres_emb_layer_;

    /// if true, solve eigenvalue problem to determine characteristic length for Nitsche's parameter
    bool nitsche_evp_;

    //! @name Fluid-fluid specific parameters
    //@{
    enum INPAR::XFEM::Monolithic_xffsi_Approach
        monolithic_approach_;  ///< type of monolithic XFFSI-approach
    enum INPAR::XFEM::XFluidFluidTimeInt xfem_timeintapproach_;  ///< XFF time integration approach
    //@}

    //! name of fluid-fluid condition for access from condition manager
    const std::string cond_name_;
  };

} /* namespace FLD */
#endif /* XFLUIDFLUID_H_ */
