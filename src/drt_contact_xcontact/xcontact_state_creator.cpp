/*----------------------------------------------------------------------------*/
/**
\file xcontact_state_creator.cpp

\brief State creator for the xcontact simulation

\maintainer Michael Hiermeier

\date Jan 19, 2017

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "xcontact_state_creator.H"
#include "xcontact_cutwizard.H"

#include "../drt_xfem/xfem_dofset.H"

#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XCONTACT::StateCreator::StateCreator()
    : XFEM::StateCreator()
{
  // empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::StateCreator::Recreate(
    Teuchos::RCP<XFEM::XFieldState> &                xstate,
    const Teuchos::RCP<DRT::DiscretizationInterface> & xfielddiscret,
    const Teuchos::RCP<DRT::DiscretizationInterface> & fielddiscret,
    const Teuchos::RCP<const Epetra_Vector> &        back_disp_col,
    const Teuchos::RCP<const Epetra_Vector> &        levelset_field_row,
    Teuchos::ParameterList &                         solver_params,
    int                                              step,
    double                                           time,
    bool                                             dosetup)
{
  // try to cast things
  Teuchos::RCP<DRT::DiscretizationXFEM> xfielddiscret_ptr =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>( xfielddiscret, true );
  Teuchos::RCP<DRT::Discretization> fielddiscret_ptr =
        Teuchos::rcp_dynamic_cast<DRT::Discretization>( fielddiscret, true );

  // if success, do the actual function call
  Recreate( xstate, xfielddiscret_ptr, fielddiscret_ptr, back_disp_col,
      levelset_field_row, solver_params, step, time, dosetup );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::StateCreator::Recreate(
    Teuchos::RCP<XFEM::XFieldState> &             xstate,
    const Teuchos::RCP<DRT::DiscretizationXFEM> & xfielddiscret,
    const Teuchos::RCP<DRT::Discretization> &     fielddiscret,
    const Teuchos::RCP<const Epetra_Vector> &     back_disp_col,
    const Teuchos::RCP<const Epetra_Vector> &     levelset_field_row,
    Teuchos::ParameterList &                      solver_params,
    int                                           step,
    double                                        time,
    bool                                          dosetup)
{
  //----------------------------------------------------------------------
  // create new cut wizard & dof-set
  Teuchos::RCP<GEO::CutWizard> wizard = Teuchos::null;
  Teuchos::RCP<XFEM::XFEMDofSet> xdofset = Teuchos::null;

   if ( CreateNewCutState( xdofset, wizard, xfielddiscret, back_disp_col,
       levelset_field_row, solver_params, step ) )
   {
     dserror( "create new cut state changed the dofset!" );
   }

  dserror("Stop, dofset did not change");
  exit( EXIT_FAILURE );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum XFEM::StateStatus XCONTACT::StateCreator::CreateNewCutState(
    Teuchos::RCP<XFEM::XFEMDofSet> &              xdofset,
    Teuchos::RCP<GEO::CutWizard> &                wizard,
    const Teuchos::RCP<DRT::DiscretizationXFEM> & xdiscret,
    const Teuchos::RCP<const Epetra_Vector> &     back_disp_col,
    const Teuchos::RCP<const Epetra_Vector> &     levelset_field_row,
    Teuchos::ParameterList &                      solver_params,
    int                                           step )
{
  CheckInitSetup();

  if ( levelset_field_row.is_null() )
    return XFEM::state_unchanged;

  //----------------------------------------------------------------------
  // create new cut wizard
  CreateCutWizard( xdiscret, back_disp_col, levelset_field_row, wizard );

  //----------------------------------------------------------------------
  // performs the "CUT"
  wizard->Cut( include_inner_ );

  //----------------------------------------------------------------------
  // set the new dofset after cut
  int max_num_my_reserved_dofs_per_node = ( maxnumdofsets_ ) * numdof_;

  // create a new XFEM-dofset
  xdofset = Teuchos::rcp( new XFEM::XFEMDofSet( *wizard ,
      max_num_my_reserved_dofs_per_node, *xdiscret ) );

  const int restart = DRT::Problem::Instance()->Restart();
  if ( (step < 1) or restart )
    minnumdofsets_ = xdiscret->DofRowMap()->MinAllGID();

  // set the minimal GID of xfem dis
  xdofset->SetMinGID( minnumdofsets_ );

 return FinishBackgroundDiscretization( xdofset, solver_params, *xdiscret );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::StateCreator::CreateCutWizard(
    const Teuchos::RCP<DRT::DiscretizationXFEM> & xdiscret,
    const Teuchos::RCP<const Epetra_Vector> &     back_disp_col,
    const Teuchos::RCP<const Epetra_Vector> &     levelset_field_row,
    Teuchos::RCP<GEO::CutWizard> &                wizard )
{
  //---------------------------------------------------------------------------
  // create a new cut wizard
  wizard = Teuchos::rcp( new XCONTACT::CutWizard( xdiscret ) );
  XCONTACT::CutWizard & xwizard = static_cast<XCONTACT::CutWizard & >( *wizard );

  //---------------------------------------------------------------------------
  // Set options for the cut wizard
  xwizard.SetOptions( INPAR::CUT::bcells_on_all_sides, nodal_dofset_strategy_,
      volume_cell_gauss_point_by_, bound_cell_gauss_point_by_, gmsh_cut_out_,
      true, false, true );

  //---------------------------------------------------------------------------
  // set background state (background mesh displacements and level-set values)
  wizard->SetBackgroundState(
      back_disp_col,
      levelset_field_row,
      1 );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum XFEM::StateStatus XCONTACT::StateCreator::FinishBackgroundDiscretization(
    const Teuchos::RCP<XFEM::XFEMDofSet> & xdofset,
    Teuchos::ParameterList &               solver_params,
    DRT::DiscretizationXFEM &              xdiscret )
{
  //---------------------------------------------------------------------------
  // check if it is necessary to create a new state object
  if ( xdiscret.IsEqualXDofSet( 0, *xdofset ) )
    return XFEM::state_unchanged;

  //field dofset has nds = 0
  xdiscret.ReplaceDofSet(0, xdofset, true );

  xdiscret.FillComplete( true, false, false );

  //print all dofsets
  xdiscret.GetDofSetProxy()->PrintAllDofsets( xdiscret.Comm() );

  //----------------------------------------------------------------------
  // recompute nullspace based on new number of dofs per node
  /* REMARK: this has to be done after replacing the discret' dofset
   * (via discret_->ReplaceDofSet) */
  xdiscret.ComputeNullSpaceIfNecessary( solver_params, true );

  return XFEM::state_changed;
}
