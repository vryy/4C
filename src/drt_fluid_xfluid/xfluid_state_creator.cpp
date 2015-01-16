/*!----------------------------------------------------------------------
\file xfluid_state_creator.cpp
\brief Creates a state object for (in)stationary XFEM fluid problems

<pre>
Maintainer:  Raffaela Kruse and Benedikt Schott
             [kruse,schott]@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15240
</pre>

*----------------------------------------------------------------------*/

#include "xfluid_state_creator.H"

#include "xfluid_state.H"
#include "xfluidfluid_state.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret_xfem.H"

#include "../linalg/linalg_mapextractor.H"

#include "../drt_cut/cut_cutwizard.H"

#include "../drt_xfem/xfem_dofset.H"
#include "../drt_xfem/xfem_neumann.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*
 |  Perform the cut and fill state container                kruse 08/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<FLD::XFluidState> FLD::XFluidStateCreator::Create(
  const Teuchos::RCP<DRT::DiscretizationXFEM>&      xdiscret,           //!< xfluid background discretization
  const Teuchos::RCP<DRT::Discretization> &         boundarydiscret,    //!< boundary discretization, a discretization whose surface elements cut the background mesh
  Teuchos::RCP<const Epetra_Vector>                 back_disp_col,      //!< col vector holding background ALE displacements for backdis
  Teuchos::RCP<const Epetra_Vector>                 cutter_disp_col,    //!< col vector holding interface displacements for cutterdis
  Teuchos::RCP<const Epetra_Vector>                 back_levelset_col,  //!< col vector holding nodal level-set values based on backdis
  std::map<int, LINALG::Matrix<3,1> > &             tip_nodes,          //!< tip nodes for crack application
  Teuchos::ParameterList &                          solver_params,
  const int                                         step,
  const double &                                    time
)
{
  //--------------------------------------------------------------------------------------
 // create new cut wizard &dofset
  Teuchos::RCP<GEO::CutWizard> wizard;
  Teuchos::RCP<XFEM::XFEMDofSet> dofset;

  CreateNewCutState(
      dofset,wizard,
      xdiscret,boundarydiscret,
      back_disp_col,cutter_disp_col,back_levelset_col,
      solver_params,
      step,
      tip_nodes);

  //--------------------------------------------------------------------------------------
  // Create the XFluidState object

  Teuchos::RCP<const Epetra_Map> xfluiddofrowmap = Teuchos::rcp(
      new Epetra_Map(*xdiscret->DofRowMap()));

  state_ = Teuchos::rcp(new FLD::XFluidState(wizard, dofset, xfluiddofrowmap));

  //--------------------------------------------------------------------------------------
  state_->SetupMapExtractors(xdiscret,time);

  return state_;
}

/*----------------------------------------------------------------------*
 |  Perform the cut and fill state container                kruse 08/14 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<FLD::XFluidFluidState> FLD::XFluidStateCreator::Create(
  const Teuchos::RCP<DRT::DiscretizationXFEM>&      xdiscret,           //!< xfluid background discretization
  const Teuchos::RCP<DRT::Discretization> &         boundarydiscret,    //!< boundary discretization, a discretization whose surface elements cut the background mesh
  const Teuchos::RCP<DRT::Discretization> &         embfluiddiscret,       //!< embedded fluid discretization
  Teuchos::RCP<const Epetra_Vector>                 back_disp_col,      //!< col vector holding background ALE displacements for backdis
  Teuchos::RCP<const Epetra_Vector>                 cutter_disp_col,    //!< col vector holding interface displacements for cutterdis
  Teuchos::RCP<const Epetra_Vector>                 back_levelset_col,  //!< col vector holding nodal level-set values based on backdis
  Teuchos::ParameterList &                          solver_params,
  const int                                         step,
  const double &                                    time
)
{
  //--------------------------------------------------------------------------------------
  // create new cut wizard & dofset
  Teuchos::RCP<GEO::CutWizard> wizard;
  Teuchos::RCP<XFEM::XFEMDofSet> dofset;

  CreateNewCutState(
      dofset,wizard,
      xdiscret,boundarydiscret,
      back_disp_col,cutter_disp_col,back_levelset_col,
      solver_params,
      step);
  //--------------------------------------------------------------------------------------
  // Create the XFluidFluidState object

  Teuchos::RCP<const Epetra_Map> xfluiddofrowmap = Teuchos::rcp(
      new Epetra_Map(*xdiscret->DofRowMap()));

  Teuchos::RCP<const Epetra_Map> embfluiddofrowmap = Teuchos::rcp(
      new Epetra_Map(*embfluiddiscret->DofRowMap()));

  Teuchos::RCP<FLD::XFluidFluidState> state = Teuchos::rcp(
      new FLD::XFluidFluidState(wizard, dofset, xfluiddofrowmap, embfluiddofrowmap));

  //--------------------------------------------------------------------------------------
  state->SetupMapExtractors(xdiscret,embfluiddiscret,time);

  state_ = state;

  return state;
}

/*----------------------------------------------------------------------*
 |  Initialize coupling matrices                           schott 12/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidStateCreator::InitCouplingMatrices(const Teuchos::RCP<DRT::Discretization>&  slavediscret)
{
  if(slavediscret == Teuchos::null) dserror("invalid slave discretization for coupling application");

  // savegraph flag set to true, as there is no change in the matrix graph expected for the lifetime
  // of this state container
  state_->C_xs_  = Teuchos::rcp(new LINALG::SparseMatrix(*state_->xfluiddofrowmap_,0,true,true,LINALG::SparseMatrix::FE_MATRIX));
  state_->C_sx_  = Teuchos::rcp(new LINALG::SparseMatrix(*slavediscret->DofRowMap(),0,true,true,LINALG::SparseMatrix::FE_MATRIX));
  state_->C_ss_  = Teuchos::rcp(new LINALG::SparseMatrix(*slavediscret->DofRowMap(),0,true,true,LINALG::SparseMatrix::FE_MATRIX));
}


/*----------------------------------------------------------------------*
 |  Initialize coupling rhs vectors                        schott 12/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidStateCreator::InitCouplingRhs(const Teuchos::RCP<DRT::Discretization>&  slavediscret)
{
  if(slavediscret == Teuchos::null) dserror("invalid slave discretization for coupling application");

  state_->rhC_s_    = LINALG::CreateVector(*slavediscret->DofRowMap(),true);
  state_->rhC_s_col_= LINALG::CreateVector(*slavediscret->DofColMap(),true);
}

/*----------------------------------------------------------------------*
 |  Initialize ALE state vectors                           schott 12/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidStateCreator::InitALEStateVectors(
    const Teuchos::RCP<DRT::DiscretizationXFEM>& xdiscret,
    Teuchos::RCP<const Epetra_Vector> dispnp_initmap,
    Teuchos::RCP<const Epetra_Vector> gridvnp_initmap
)
{

  //! @name Ale Displacement at time n+1
  state_->dispnp_ = LINALG::CreateVector(*state_->xfluiddofrowmap_,true);
  xdiscret->ExportInitialtoActiveVector(dispnp_initmap,state_->dispnp_);

  //! @name Grid Velocity at time n+1
  state_->gridvnp_ = LINALG::CreateVector(*state_->xfluiddofrowmap_,true);
  xdiscret->ExportInitialtoActiveVector(gridvnp_initmap,state_->gridvnp_);
}

/*----------------------------------------------------------------------*
 |  Initialize ALE state vectors                           schott 12/14 |
 *----------------------------------------------------------------------*/
void FLD::XFluidStateCreator::CreateNewCutState(
  Teuchos::RCP<XFEM::XFEMDofSet> &                  dofset,             //!< xfem dofset obtained from the new wizard
  Teuchos::RCP<GEO::CutWizard> &                    wizard,             //!< cut wizard associated with current intersection state
  const Teuchos::RCP<DRT::DiscretizationXFEM>&      xdiscret,           //!< xfluid background discretization
  const Teuchos::RCP<DRT::Discretization> &         boundarydiscret,    //!< boundary discretization, a discretization whose surface elements cut the background mesh
  Teuchos::RCP<const Epetra_Vector>                 back_disp_col,      //!< col vector holding background ALE displacements for backdis
  Teuchos::RCP<const Epetra_Vector>                 cutter_disp_col,    //!< col vector holding interface displacements for cutterdis
  Teuchos::RCP<const Epetra_Vector>                 back_levelset_col,  //!< col vector holding nodal level-set values based on backdis
  Teuchos::ParameterList &                          solver_params,
  const int                                         step,
  std::map<int, LINALG::Matrix<3,1> >               tip_nodes           //!< tip nodes for crack application
)
{
  // new wizard
  wizard = Teuchos::rcp( new GEO::CutWizard(xdiscret, boundarydiscret) );

  // Set options for the cut wizard
  wizard->SetOptions(
      VolumeCellGaussPointBy_,       // how to create volume cell Gauss points?
      BoundCellGaussPointBy_,        // how to create boundary cell Gauss points?
      gmsh_cut_out_,                 // gmsh output for cut library
      true,                          // find point positions
      false,                         // generate only tet cells
      true                           // print screen output
  );

  // Set the state vectors used for the cut
  wizard->SetState(
      back_disp_col,      //!< col vector holding background ALE displacements for backdis
      cutter_disp_col,    //!< col vector holding interface displacements for cutterdis
      back_levelset_col   //!< col vector holding nodal level-set values based on backdis
  );

  // Set crack tip nodes in case of crack application
  wizard->setCrackTipNodes(tip_nodes);

  //--------------------------------------------------------------------------------------
  // performs the "CUT"
  wizard->Cut(include_inner_);
  //--------------------------------------------------------------------------------------


  //--------------------------------------------------------------------------------------
  // set the new dofset after cut
  int maxNumMyReservedDofsperNode = (maxnumdofsets_)*4;

  // create a new XFEM-dofset
  dofset = Teuchos::rcp( new XFEM::XFEMDofSet( wizard , maxNumMyReservedDofsperNode, xdiscret ) );

  const int restart = DRT::Problem::Instance()->Restart();
  if ((step < 1) or restart)
    minnumdofsets_ = xdiscret->DofRowMap()->MinAllGID();

  dofset->SetMinGID(minnumdofsets_); // set the minimal GID of xfem dis
  xdiscret->ReplaceDofSet( dofset, true );

  xdiscret->FillComplete(true,false,false);

  //print all dofsets
  xdiscret->GetDofSetProxy()->PrintAllDofsets(xdiscret->Comm());

  //--------------------------------------------------------------------------------------
  // recompute nullspace based on new number of dofs per node
  // REMARK: this has to be done after replacing the discret' dofset (via discret_->ReplaceDofSet)
  xdiscret->ComputeNullSpaceIfNecessary(solver_params,true);
}
