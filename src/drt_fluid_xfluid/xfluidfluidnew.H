/*!----------------------------------------------------------------------
\file xfluidfluidnew.H
\brief Control routine for fluid-fluid (in)stationary solvers with XFEM

ATTENTION! Class is still a prototype. Does not provide full (and correct)
functionality yet!

<pre>
Maintainer:  Raffaela Kruse
             kruse@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>
*----------------------------------------------------------------------*/

#ifndef XFLUIDFLUIDNEW_H_
#define XFLUIDFLUIDNEW_H_

#include "xfluid.H"
#include "xfluidfluid_state.H"

#include "../drt_fluid/fluidimplicitintegration.H"

namespace XFEM
{
  class XFluidFluidTimeIntegration;
}

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
class XFluidFluidNew : public XFluid
{
public:

  //! ctor
  XFluidFluidNew(
    const Teuchos::RCP<FLD::FluidImplicitTimeInt> & embedded_fluid,     ///< embedded fluid
    const Teuchos::RCP<DRT::Discretization>&        xfluiddis,          ///< background fluid discretization
    const Teuchos::RCP<LINALG::Solver>&             solver,             ///< fluid solver
    const Teuchos::RCP<Teuchos::ParameterList>&     params,             ///< xfluid params
    const Teuchos::RCP<IO::DiscretizationWriter>&   output,             ///< discretization writer for paraview output
    bool                                            ale_xfluid = false, ///< background (XFEM) fluid in ALE-formulation
    bool                                            ale_fluid  = false  ///< embedded fluid in ALE-formulation
  );

  virtual ~XFluidFluidNew();

  /// initialization
  virtual void Init();

  /// set fluid-fluid coupling specific parameters
  void SetXFluidFluidParams();

  /// set initial flow field for fluid domains
  void SetInitialFlowField(
    const INPAR::FLUID::InitialField initfield,
    const int startfuncno
  )
  {
    XFluid::SetInitialFlowField(initfield,startfuncno);
    embedded_fluid_->SetInitialFlowField(initfield,startfuncno);
  }

  /// Update the solution after convergence of the nonlinear
  /// iteration. Current solution becomes old solution of next timestep.
  void TimeUpdate();

  /// CUT at new interface position, transform vectors,
  /// perform time integration and set new Vectors
  void CutAndSetStateVectors();

  /*!
   * \brief set underlying dof-maps for new shape derivatives matrix
   * (required in the context of monolithic XFFSI)
   * \param (in) fsiextractor map extractor to distinguish between FSI/non-FSI dof
   * in the merged fluid-fluid-dof-map
   * \param (in) condelements map of conditioned elements
   */
  void PrepareShapeDerivatives(
      const Teuchos::RCP<const LINALG::MultiMapExtractor> fsiextractor,
      const Teuchos::RCP<std::set<int> > condelements)
  {
    dserror("Not implemented yet!");
  }

  /*!
   * \brief return fluid-fluid system matrix as block matrix
   * (required in the context of monolithic fluidsplit-XFFSI)
   * \param (in) innermap map of inner embedded and background fluid dof
   * \param (in) condmap map of fsi interface dof
   * \return coupled fluid-fluid block system matrix
   */
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix(
      Teuchos::RCP<Epetra_Map> innermap,
      Teuchos::RCP<Epetra_Map> condmap);

  Teuchos::RCP<const Epetra_Vector> InitialGuess()  { return state_->xffluidincvel_; }
  Teuchos::RCP<const Epetra_Vector> RHS()           { return state_->xffluidresidual_; }

  Teuchos::RCP<const Epetra_Vector> Velnp()             { return state_->xffluidvelnp_; }
  Teuchos::RCP<      Epetra_Vector> WriteAccessVelnp()  { return state_->xffluidvelnp_; }
  Teuchos::RCP<const Epetra_Vector> Veln()              { return state_->xffluidveln_; }
  Teuchos::RCP<const Epetra_Vector> GridVel()           { return embedded_fluid_->GridVel(); }
  Teuchos::RCP<      Epetra_Vector> WriteAccessGridVel(){ return embedded_fluid_->WriteAccessGridVel(); }
  Teuchos::RCP<const Epetra_Vector> Stepinc()           { return stepinc_; }

  Teuchos::RCP<const Epetra_Vector> Dispnp()                  { return embedded_fluid_->Dispnp(); }
  Teuchos::RCP<      Epetra_Vector> WriteAccessDispnp()       { return embedded_fluid_->WriteAccessDispnp(); }
  Teuchos::RCP<      Epetra_Vector> WriteAccessDispOldState() { return dispnpoldstate_; }
  Teuchos::RCP<const Epetra_Vector> Dispn()                   { return embedded_fluid_->Dispn(); }

  Teuchos::RCP<const Epetra_Map> DofRowMap()
  {
    return state_->xffluiddofrowmap_;
  }

  Teuchos::RCP<const Epetra_Map> PressureRowMap();
  Teuchos::RCP<const Epetra_Map> VelocityRowMap();

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives(){ return extended_shapederivatives_; }

  /// get the combined fluid-fluid system matrix
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix()
  {
    return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(state_->xffluidsysmat_);
  }

  /// get map extractor for background/embedded fluid

  Teuchos::RCP<FLD::UTILS::XFluidFluidMapExtractor>const& XFluidFluidMapExtractor(){ return state_->xffluidsplitter_; }

  //! get combined background and embedded fluid dirichlet map extractor
  Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() { return state_->xffluiddbcmaps_; }

  //! access to new merged fluid state container
  Teuchos::RCP<FLD::XFluidFluidState> GetNewState();

private:

  ///  prepare embedded discretization for coupling type embedded-sided Nitsche
  void PrepareEmbeddedDistribution();

  ///  prepare mapping between boundary elements and outer embedded elements
  /// (required for coupling type embedded-sided Nitsche)
  void CreateBoundaryEmbeddedMap();

  /// create new discretization from all embedded fluid elements,
  /// that contribute to the fluid-fluid interface with at least one face
  void CreateEmbeddedBoundaryDiscretization();

  //! create a new merged state container
  void CreateState();

  /// call the loop over elements to assemble volume and interface integrals
  void AssembleMatAndRHS(
      int itnum                           ///< iteration number
  );

  /// embedded fluid field
  Teuchos::RCP<FLD::FluidImplicitTimeInt> embedded_fluid_;

  /// whether the embedded fluid is an ALE-fluid
  bool alefluid_;

  /// auxiliary discretization of outer embedded elements
  Teuchos::RCP<DRT::Discretization>       embboundarydis_;

  /// merged fluid state (cut-dependent) at time n+1
  Teuchos::RCP<FLD::XFluidFluidState> state_;
  /// merged fluid state (cut-dependent) at time n
  Teuchos::RCP<FLD::XFluidFluidState> staten_;

  /// vector with Newton increments (used for monolithic fully newton fsi approach)
  Teuchos::RCP<Epetra_Vector>         stepinc_;

  /// ALE-displacements of previous step
  Teuchos::RCP<Epetra_Vector>         dispnpoldstate_;

  /// shape derivatives matrix (linearization with respect to mesh motion),
  /// including background fluid dof (set to zero)
  Teuchos::RCP<LINALG::BlockSparseMatrixBase>   extended_shapederivatives_;

  /// xfem time integration object
  Teuchos::RCP<XFEM::XFluidFluidTimeIntegration> xfluidfluid_timeint_;

  /// flag, that indicates active shape derivatives
  bool active_shapederivatives_;

  /// flag, that indicates evaluation of edge-based stabilizing terms
  bool eval_eos_;

  /// flag, that indicates evaluation of pressure EOS-terms on outer embedded elements
  bool xff_eos_pres_emb_layer_;

  /// if true, solve eigenvalue problem to determine characteristic length for Nitsche's parameter
  bool nitsche_evp_;

  /// for coupling using stress-based LM: L2-projection between stress fields on whole cut element
  /// or on physical volume
  enum INPAR::XFEM::Hybrid_LM_L2_Proj hybrid_lm_l2_proj_;

  /// convective stabilization scaling for inflow terms
  enum INPAR::XFEM::XFF_ConvStabScaling xff_conv_stab_scaling_;

  //! @name Fluid-fluid specific parameters
  //@{
  enum INPAR::XFEM::Monolithic_xffsi_Approach monolithic_approach_;  ///< type of monolithic XFFSI-approach
  enum INPAR::XFEM::XFluidFluidTimeInt        xfem_timeintapproach_; ///< XFF time integration approach
  //@}

};

} /* namespace FLD */
#endif /* XFLUIDFLUIDNEW_H_ */
