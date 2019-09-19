/*----------------------------------------------------------------------*/
/*! \file

\brief State class for (in)stationary XFEM fluid problems

\level 0

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
 */
/*----------------------------------------------------------------------*/

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_utils.H"

#include "../drt_cut/cut_cutwizard.H"

#include "../drt_xfem/xfem_condition_manager.H"
#include "../drt_xfem/xfem_dofset.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_xfem/xfield_state_utils.H"
#include "../drt_fluid/fluid_utils.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "xfluid_state.H"


/*----------------------------------------------------------------------*
 |  ctor  Initialize coupling matrices                     schott 01/15 |
 *----------------------------------------------------------------------*/
FLD::XFluidState::CouplingState::CouplingState(
    const Teuchos::RCP<const Epetra_Map>& xfluiddofrowmap,
    const Teuchos::RCP<DRT::Discretization>& slavediscret_mat,
    const Teuchos::RCP<DRT::Discretization>& slavediscret_rhs)
    : is_active_(true)
{
  if (slavediscret_mat == Teuchos::null)
    dserror("invalid slave discretization for coupling application");
  if (slavediscret_rhs == Teuchos::null)
    dserror("invalid slave discretization for coupling application");


  // savegraph flag set to true, as there is no change in the matrix graph expected for the lifetime
  // of this state container
  // no explicit Dirichlet, otherwise new matrices will be created in ApplyDirichlet
  // NOTE: setting explicit Dirichlet to false can cause problems with ML preconditioner (see remark
  // in LINALG::Sparsematrix) however, we prefer not to build new matrices in ApplyDirichlet
  C_xs_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *xfluiddofrowmap, 300, false, true, LINALG::SparseMatrix::FE_MATRIX));
  C_sx_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *slavediscret_mat->DofRowMap(), 300, false, true, LINALG::SparseMatrix::FE_MATRIX));
  C_ss_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *slavediscret_mat->DofRowMap(), 300, false, true, LINALG::SparseMatrix::FE_MATRIX));

  rhC_s_ = LINALG::CreateVector(*slavediscret_rhs->DofRowMap(), true);
  rhC_s_col_ = LINALG::CreateVector(*slavediscret_rhs->DofColMap(), true);
}

/*----------------------------------------------------------------------*
 |  zero coupling matrices and rhs vectors                 schott 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::CouplingState::ZeroCouplingMatricesAndRhs()
{
  if (!is_active_) return;

  // zero all coupling matrices and rhs vectors
  XFEM::ZeroMatrix(C_xs_);
  XFEM::ZeroMatrix(C_sx_);
  XFEM::ZeroMatrix(C_ss_);

  rhC_s_->PutScalar(0.0);
  rhC_s_col_->PutScalar(0.0);
}

/*----------------------------------------------------------------------*
 |  complete coupling matrices and rhs vectors             schott 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::CouplingState::CompleteCouplingMatricesAndRhs(
    const Epetra_Map& xfluiddofrowmap, const Epetra_Map& slavedofrowmap)
{
  if (!is_active_) return;

  // REMARK: for EpetraFECrs matrices Complete() calls the GlobalAssemble() routine to gather
  // entries from all processors (domain-map are the columns, range-map are the rows)
  C_xs_->Complete(slavedofrowmap, xfluiddofrowmap);
  C_sx_->Complete(xfluiddofrowmap, slavedofrowmap);
  C_ss_->Complete(slavedofrowmap, slavedofrowmap);

  //  std::cout << "number of nonzeros: C_xs" << C_xs_->EpetraMatrix()->MaxNumEntries() <<
  //  std::endl; std::cout << "number of nonzeros: C_sx" << C_sx_->EpetraMatrix()->MaxNumEntries()
  //  << std::endl; std::cout << "number of nonzeros: C_ss" <<
  //  C_ss_->EpetraMatrix()->MaxNumEntries() << std::endl;
  //-------------------------------------------------------------------------------
  // export the rhs coupling vector to a row vector
  Epetra_Vector rhC_s_tmp(rhC_s_->Map(), true);
  Epetra_Export exporter_rhC_s_col(rhC_s_col_->Map(), rhC_s_tmp.Map());
  int err = rhC_s_tmp.Export(*rhC_s_col_, exporter_rhC_s_col, Add);
  if (err) dserror("Export using exporter returned err=%d", err);

  rhC_s_->Update(1.0, rhC_s_tmp, 0.0);
}


/*----------------------------------------------------------------------*
 |  destroy the coupling objects and it's content          schott 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::CouplingState::Destroy(bool throw_exception)
{
  if (!is_active_) return;

  XFEM::DestroyMatrix(C_xs_, throw_exception);
  XFEM::DestroyMatrix(C_sx_, throw_exception);
  XFEM::DestroyMatrix(C_ss_, throw_exception);

  XFEM::DestroyRCPObject(rhC_s_, throw_exception);
  XFEM::DestroyRCPObject(rhC_s_col_, throw_exception);

  is_active_ = false;
}


/*----------------------------------------------------------------------*
 |  Constructor for XFluidState                             kruse 08/14 |
 *----------------------------------------------------------------------*/
FLD::XFluidState::XFluidState(const Teuchos::RCP<XFEM::ConditionManager>& condition_manager,
    const Teuchos::RCP<GEO::CutWizard>& wizard, const Teuchos::RCP<XFEM::XFEMDofSet>& dofset,
    const Teuchos::RCP<const Epetra_Map>& xfluiddofrowmap,
    const Teuchos::RCP<const Epetra_Map>& xfluiddofcolmap)
    : xfluiddofrowmap_(xfluiddofrowmap),
      xfluiddofcolmap_(xfluiddofcolmap),
      dofset_(dofset),
      wizard_(wizard),
      condition_manager_(condition_manager)
{
  InitSystemMatrix();
  InitStateVectors();

  InitCouplingMatricesAndRhs();
}

/*----------------------------------------------------------------------*
 |  Initialize (coupling) matrices & rhs-vectors
 |                                                         kruse 08/14
 *----------------------------------------------------------------------*/
void FLD::XFluidState::InitSystemMatrix()
{
  // create an EpetraFECrs matrix that does communication for non-local rows and columns
  // * this enables to do the evaluate loop over just row elements instead of col elements
  // * time consuming assemble for cut elements is done only once on a unique row processor
  // REMARK: call the SparseMatrix: * explicitdirichlet = false (is used in ApplyDirichlet, true
  // would create a new matrix when DBS will be applied)
  //                                    setting flag to false can cause problems with ML
  //                                    preconditioner (see remark in LINALG::Sparsematrix) however,
  //                                    we prefer not to build new matrices in ApplyDirichlet
  //                                * savegraph = true/false: To save the graph (pattern for
  //                                non-zero entries) leads to a speedup in the assembly of the
  //                                matrix
  //                                    for subsequent assemblies.
  //                                    However, do not use the savegraph-option when the matrix
  //                                    graph can change during the usage of this object of
  //                                    SparseMatrix or use the Reset()-function instead of Zero().
  //                                    For XFEM-problems, the matrix graph changes between
  //                                    timesteps, however, then a new state class and sparsematrix
  //                                    is created, otherwise a Reset()-function has to be called
  //                                    instead of the Zero()-function. We are using the save-graph
  //                                    option.
  // * the estimate of the number of nonzero entries is adapted to hex8 elements with 8 adjacent
  // elements around a node
  //   + edge-based couplings component-wise v_x->u_x, v_y->u_y, v_z->u_z, q->p
  //   number of non-zeros (for hex8 elements): 108+54 = 162
  sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *xfluiddofrowmap_, 162, false, true, LINALG::SparseMatrix::FE_MATRIX));
}


/*----------------------------------------------------------------------*
 |  Initialize state vectors                                kruse 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::InitStateVectors()
{
  // Vectors passed to the element
  // -----------------------------
  // velocity/pressure at time n+1, n and n-1
  velnp_ = LINALG::CreateVector(*xfluiddofrowmap_, true);
  veln_ = LINALG::CreateVector(*xfluiddofrowmap_, true);
  velnm_ = LINALG::CreateVector(*xfluiddofrowmap_, true);

  // velocity/pressure at time n+alpha_F
  velaf_ = LINALG::CreateVector(*xfluiddofrowmap_, true);


  // acceleration/(scalar time derivative) at time n+1 and n
  accnp_ = LINALG::CreateVector(*xfluiddofrowmap_, true);
  accn_ = LINALG::CreateVector(*xfluiddofrowmap_, true);

  // acceleration/(scalar time derivative) at time n+alpha_M/(n+alpha_M/n)
  accam_ = LINALG::CreateVector(*xfluiddofrowmap_, true);

  // scalar at time n+alpha_F/n+1 and n+alpha_M/n
  // (only required for low-Mach-number case)
  // ... this is a dummy to avoid errors
  scaaf_ = LINALG::CreateVector(*xfluiddofrowmap_, true);
  scaam_ = LINALG::CreateVector(*xfluiddofrowmap_, true);

  // history vector
  hist_ = LINALG::CreateVector(*xfluiddofrowmap_, true);

  // the vector containing body and surface forces
  neumann_loads_ = LINALG::CreateVector(*xfluiddofrowmap_, true);

  // rhs: standard (stabilized) residual vector (rhs for the incremental form)
  residual_ = LINALG::CreateVector(*xfluiddofrowmap_, true);
  residual_col_ = LINALG::CreateVector(*xfluiddofcolmap_, true);
  trueresidual_ = LINALG::CreateVector(*xfluiddofrowmap_, true);

  // nonlinear iteration increment vector
  incvel_ = LINALG::CreateVector(*xfluiddofrowmap_, true);

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_ = LINALG::CreateVector(*xfluiddofrowmap_, true);
}

/*----------------------------------------------------------------------*
 |  Initialize coupling matrices                           schott 01/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::InitCouplingMatricesAndRhs()
{
  // loop all coupling objects
  for (int coup_idx = 0; coup_idx < condition_manager_->NumCoupling(); coup_idx++)
  {
    Teuchos::RCP<XFEM::CouplingBase> coupling = condition_manager_->GetCouplingByIdx(coup_idx);

#ifdef DEBUG
    if (coupling == Teuchos::null) dserror("invalid coupling object!");
#endif

    Teuchos::RCP<XFluidState::CouplingState> coup_state = Teuchos::null;

    if (!condition_manager_->IsCouplingCondition(
            coupling->GetName()))  // coupling or one-sided non-coupling object
    {
      // create coupling state object with coupling matrices initialized with Teuchos::null
      coup_state = Teuchos::rcp(new XFluidState::CouplingState());
    }
    else
    {
      if (condition_manager_->IsLevelSetCondition(coup_idx))
      {
        // coupling matrices can be assembled into the fluid sysmat
        // coupling rhs terms can be assembled into the fluid residual

        // create weak rcp's which simplifies deleting the systemmatrix
        Teuchos::RCP<LINALG::SparseMatrix> sysmat_weakRCP =
            sysmat_.create_weak();  // no increment in strong reference counter
        Teuchos::RCP<Epetra_Vector> residual_weakRCP = residual_.create_weak();
        Teuchos::RCP<Epetra_Vector> residual_col_weakRCP = residual_col_.create_weak();

        coup_state = Teuchos::rcp(new XFluidState::CouplingState(sysmat_weakRCP, sysmat_weakRCP,
            sysmat_weakRCP, residual_weakRCP, residual_col_weakRCP));
      }
      else if (condition_manager_->IsMeshCondition(coup_idx))
      {
        // for matrix use the full condition dis to enable assign in blockmatrix (row map of matrix
        // is not changed by complete!) for rhs we need the additional ghosting of the boundary zone
        // on slave side, therefore the specifically modified coupling dis
        coup_state = Teuchos::rcp(new XFluidState::CouplingState(
            xfluiddofrowmap_, coupling->GetCondDis(), coupling->GetCouplingDis()));
      }
      else
        dserror(
            "coupling object is neither a level-set coupling object nor a mesh-coupling object");
    }

    coup_state_[coup_idx] = coup_state;

  }  // loop coupling objects
}

/*----------------------------------------------------------------------*
 |  Initialize ALE state vectors                           schott 12/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::InitALEStateVectors(const Teuchos::RCP<DRT::DiscretizationXFEM>& xdiscret,
    Teuchos::RCP<const Epetra_Vector> dispnp_initmap,
    Teuchos::RCP<const Epetra_Vector> gridvnp_initmap)
{
  //! @name Ale Displacement at time n+1
  dispnp_ = LINALG::CreateVector(*xfluiddofrowmap_, true);
  xdiscret->ExportInitialtoActiveVector(dispnp_initmap, dispnp_);

  //! @name Grid Velocity at time n+1
  gridvnp_ = LINALG::CreateVector(*xfluiddofrowmap_, true);
  xdiscret->ExportInitialtoActiveVector(gridvnp_initmap, gridvnp_);
}


/*----------------------------------------------------------------------*
 |  Initialize state vectors                               schott 12/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::ZeroCouplingMatricesAndRhs()
{
  // loop all coupling objects
  for (std::map<int, Teuchos::RCP<CouplingState>>::iterator cs = coup_state_.begin();
       cs != coup_state_.end(); cs++)
    cs->second->ZeroCouplingMatricesAndRhs();
}

/*----------------------------------------------------------------------*
 |  Complete coupling matrices and rhs vectors             schott 12/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::CompleteCouplingMatricesAndRhs()
{
  CompleteCouplingMatricesAndRhs(xfluiddofrowmap_);
}


/*----------------------------------------------------------------------*
 |  Complete coupling matrices and rhs vectors             schott 12/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::CompleteCouplingMatricesAndRhs(
    const Teuchos::RCP<const Epetra_Map>& fluiddofrowmap  ///< fluid dof row map used for complete
)
{
  // loop all coupling objects
  for (std::map<int, Teuchos::RCP<CouplingState>>::iterator cs = coup_state_.begin();
       cs != coup_state_.end(); cs++)
  {
    int coupl_idx = cs->first;

    // complete only the mesh coupling objects, levelset-couplings are completed by the
    // sysmat->Complete call
    if (condition_manager_->IsMeshCondition(coupl_idx))
    {
      Teuchos::RCP<XFEM::CouplingBase> coupling = condition_manager_->GetCouplingByIdx(coupl_idx);
      cs->second->CompleteCouplingMatricesAndRhs(
          *fluiddofrowmap, *coupling->GetCondDis()->DofRowMap());  // complete w.r.t
    }
  }
}


/*----------------------------------------------------------------------*
 |    /// zero system matrix and related rhs vectors       schott 08/15 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::ZeroSystemMatrixAndRhs()
{
  XFEM::ZeroMatrix(sysmat_);

  // zero residual vectors
  residual_col_->PutScalar(0.0);
  residual_->PutScalar(0.0);
  trueresidual_->PutScalar(0.0);
}


/*----------------------------------------------------------------------*
 |  Set dirichlet- and velocity/pressure-map extractor      kruse 08/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidState::SetupMapExtractors(
    const Teuchos::RCP<DRT::Discretization>& xfluiddiscret, const double& time)
{
  // create dirichlet map extractor
  Teuchos::ParameterList eleparams;
  // other parameters needed by the elements
  eleparams.set("total time", time);
  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  xfluiddiscret->EvaluateDirichlet(
      eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);

  zeros_->PutScalar(0.0);

  // create vel-pres splitter
  const int numdim = DRT::Problem::Instance()->NDim();
  velpressplitter_ = Teuchos::rcp(new LINALG::MapExtractor());
  FLD::UTILS::SetupFluidSplit(*xfluiddiscret, numdim, 1, *velpressplitter_);
}


bool FLD::XFluidState::Destroy()
{
  // destroy system matrix (destroy after coupling matrices, as for twophase problems the coupling
  // matrices are identical to the system matrix)
  XFEM::DestroyMatrix(sysmat_);

  // destroy all coupling system matrices and rhs vectors (except for levelset coupling objects
  for (std::map<int, Teuchos::RCP<CouplingState>>::iterator i = coup_state_.begin();
       i != coup_state_.end(); ++i)
  {
    Teuchos::RCP<CouplingState>& cs = i->second;  // RCP reference!
    if (cs != Teuchos::null)
    {
      cs->Destroy(true);   // throw exception when object could not be deleted or more than one
                           // strong RCP points to it
      cs = Teuchos::null;  // invalidate coupling state object
    }
  }

  // destroy dofrowmap and dofcolmap
  XFEM::DestroyRCPObject(xfluiddofrowmap_);
  XFEM::DestroyRCPObject(xfluiddofcolmap_);

  // destroy state vectors
  XFEM::DestroyRCPObject(velnp_);
  XFEM::DestroyRCPObject(veln_);
  XFEM::DestroyRCPObject(velnm_);
  XFEM::DestroyRCPObject(velaf_);

  XFEM::DestroyRCPObject(accnp_);
  XFEM::DestroyRCPObject(accn_);
  XFEM::DestroyRCPObject(accam_);

  XFEM::DestroyRCPObject(scaaf_);
  XFEM::DestroyRCPObject(scaam_);

  XFEM::DestroyRCPObject(hist_);
  XFEM::DestroyRCPObject(neumann_loads_);

  XFEM::DestroyRCPObject(residual_);
  XFEM::DestroyRCPObject(trueresidual_);

  XFEM::DestroyRCPObject(zeros_);
  XFEM::DestroyRCPObject(incvel_);

  XFEM::DestroyRCPObject(dispnp_);
  XFEM::DestroyRCPObject(gridvnp_);


  // destroy velpressplitter_
  XFEM::DestroyRCPObject(velpressplitter_);

  // wizard, dofset and conditionmanager keep RCPs pointing to them and cannot be destroyed as they
  // are further used in xfluid-class decrease at least the strong reference counter
  wizard_ = Teuchos::null;
  dofset_ = Teuchos::null;
  condition_manager_ = Teuchos::null;

  return true;
}

void FLD::XFluidState::UpdateBoundaryCellCoords()
{
  // loop all mesh coupling objects
  for (int mc_idx = 0; mc_idx < condition_manager_->NumMeshCoupling(); mc_idx++)
  {
    Teuchos::RCP<XFEM::MeshCoupling> mc_coupl = condition_manager_->GetMeshCoupling(mc_idx);

    if (!mc_coupl->CutGeometry()) continue;  // If don't cut the background mesh.

    wizard_->UpdateBoundaryCellCoords(mc_coupl->GetCutterDis(), mc_coupl->GetCutterDispCol(),
        condition_manager_->GetMeshCouplingStartGID(mc_idx));
  }
}
