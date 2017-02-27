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

#include "../drt_adapter/ad_str_structure_new.H"

#include "../drt_xfem/xfem_dofset.H"
#include "../drt_xfem/xfield_state.H"

#include "../drt_lib/drt_discret_xfem.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io_pstream.H"
#include <Teuchos_Time.hpp>

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
    const Teuchos::RCP<ADAPTER::Field> &             xfield,
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
  Recreate( xstate, xfield, xfielddiscret_ptr, fielddiscret_ptr, back_disp_col,
      levelset_field_row, solver_params, step, time, dosetup );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::StateCreator::Recreate(
    Teuchos::RCP<XFEM::XFieldState> &             xstate,
    const Teuchos::RCP<ADAPTER::Field> &          xfield,
    const Teuchos::RCP<DRT::DiscretizationXFEM> & xfielddiscret,
    const Teuchos::RCP<DRT::Discretization> &     fielddiscret,
    const Teuchos::RCP<const Epetra_Vector> &     back_disp_col,
    const Teuchos::RCP<const Epetra_Vector> &     levelset_field_row,
    Teuchos::ParameterList &                      solver_params,
    int                                           step,
    double                                        time,
    bool                                          dosetup)
{
  //---------------------------------------------------------------------------
  // create new cut wizard & dof-set
  Teuchos::RCP<GEO::CutWizard> wizard = Teuchos::null;
  Teuchos::RCP<XFEM::XFEMDofSet> xdofset = Teuchos::null;

  enum XFEM::StateStatus status = CreateNewCutState( xdofset, wizard,
      xfielddiscret, back_disp_col, levelset_field_row, solver_params, step );

  IO::cout << "\nXCONTACT::StateCreator::Recreate:" << IO::endl;

  switch ( status )
  {
    case XFEM::state_changed:
    {
      Teuchos::RCP<XFEM::XFieldState> xstate_new = Teuchos::null;
      {
        const double t_start = Teuchos::Time::wallTime();
        IO::cout << "\t* 1/3 CreateNewXFieldState ...";

        xstate_new = CreateNewXFieldState( xfield,
            wizard, xdofset, xfielddiscret, fielddiscret );

        const double t_diff = Teuchos::Time::wallTime() - t_start;
        IO::cout << " Success (" << t_diff << " secs)" << IO::endl;
      }

      {
        const double t_start = Teuchos::Time::wallTime();
        IO::cout << "\t* 2/3 DestroyXFieldState .....";

        DestroyXFieldState( xfield );

        const double t_diff = Teuchos::Time::wallTime() - t_start;
        IO::cout << " Success (" << t_diff << " secs)" << IO::endl;
      }

      {
        const double t_start = Teuchos::Time::wallTime();
        IO::cout << "\t* 3/3 SetNewState ............";

        xstate->SetNewState( *xstate_new );

        const double t_diff = Teuchos::Time::wallTime() - t_start;
        IO::cout << " Success (" << t_diff << " secs)" << IO::endl;
      }

      break;
    }
    case XFEM::state_unchanged:
    {
      {
        const double t_start = Teuchos::Time::wallTime();
        IO::cout << "\t* 1/1 ResetXFieldNonStandardDofs ...";

        ResetXFieldNonStandardDofs( xfield );

        const double t_diff = Teuchos::Time::wallTime() - t_start;
        IO::cout << " Success (" << t_diff << " secs)" << IO::endl;
      }

      break;
    }
    default:
    {
      dserror( "Unsupported StateStatus: %s",
          XFEM::StateStatus2String( status ).c_str() );
      exit( EXIT_FAILURE );
    }
  }
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

  //---------------------------------------------------------------------------
  // create new cut wizard
  CreateCutWizard( xdiscret, back_disp_col, levelset_field_row, wizard );

  //---------------------------------------------------------------------------
  // performs the "CUT"
  wizard->Cut( include_inner_ );

  //---------------------------------------------------------------------------
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

  // field dofset has nds = 0
  xdiscret.ReplaceDofSet(0, xdofset, true );

  xdiscret.FillComplete( true, false, false );

  //print all dofsets
  xdiscret.GetDofSetProxy()->PrintAllDofsets( xdiscret.Comm() );

  //---------------------------------------------------------------------------
  // recompute nullspace based on new number of dofs per node
  /* REMARK: this has to be done after replacing the discret' dofset
   * (via discret_->ReplaceDofSet) */
  xdiscret.ComputeNullSpaceIfNecessary( solver_params, true );

  return XFEM::state_changed;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<XFEM::XFieldState> XCONTACT::StateCreator::CreateNewXFieldState(
    const Teuchos::RCP<ADAPTER::Field> &   xfield,
    const Teuchos::RCP<GEO::CutWizard> &   wizard,
    const Teuchos::RCP<XFEM::XFEMDofSet> & xdofset,
    const Teuchos::RCP<DRT::DiscretizationInterface> & xfielddiscret,
    const Teuchos::RCP<DRT::DiscretizationInterface> & fielddiscret ) const
{
  Teuchos::RCP<XFEM::XFieldState> new_xstate = CreateXFieldFieldState();
  new_xstate->Init( Teuchos::null, wizard, xdofset, xfielddiscret, fielddiscret );

  Teuchos::RCP<ADAPTER::StructureNew> xstructure =
      Teuchos::rcp_dynamic_cast<ADAPTER::StructureNew>( xfield, true );

  xstructure->CreateNewXFieldState( new_xstate );

  return new_xstate;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::StateCreator::DestroyXFieldState(
    const Teuchos::RCP<ADAPTER::Field> & xfield ) const
{
  Teuchos::RCP<ADAPTER::StructureNew> xstructure =
      Teuchos::rcp_dynamic_cast<ADAPTER::StructureNew>( xfield, true );

  xstructure->DestroyState();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XCONTACT::StateCreator::ResetXFieldNonStandardDofs(
    const Teuchos::RCP<ADAPTER::Field> & xfield ) const
{
  Teuchos::RCP<ADAPTER::StructureNew> xstructure =
      Teuchos::rcp_dynamic_cast<ADAPTER::StructureNew>( xfield, true );

  xstructure->ResetXFieldNonStandardDofs();
}
