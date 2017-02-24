/*----------------------------------------------------------------------------*/
/**
\file xstr_xstructure_state.cpp

\brief XStructure state handling ( only one single XFEM-discretization )

\maintainer Michael Hiermeier

\date Jun 28, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "xstr_xstructure_state.H"

#include "../drt_structure_new/str_utils.H"

#include "../drt_xfem/xfield_state_utils.H"

#include "../drt_lib/drt_discret_interface.H"
#include "../drt_lib/drt_node.H"

#include <Epetra_Export.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XSTR::XStructureState::XStructureState()
{
  // intentionally left blank
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::XStructureState::Setup()
{
  CheckInit();

  STR::TIMINT::BaseDataGlobalState::Setup();

  RegisterStateVectorsInMap();
  RegisterStateMatricesInMap();

  XFEM::XFieldState::issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::XStructureState::RegisterStateVectorsInMap()
{
  xstate_vectors_[ state_vec_disnp ] = & GetMutableDisNp();
  xstate_vectors_[ state_vec_velnp ] = & GetMutableVelNp();
  xstate_vectors_[ state_vec_accnp ] = & GetMutableAccNp();

  xstate_vectors_[ state_vec_fintn ] = & GetMutableFintN();
  xstate_vectors_[ state_vec_fintnp ] = & GetMutableFintNp();
  xstate_vectors_[ state_vec_fextn ] = & GetMutableFextN();
  xstate_vectors_[ state_vec_fextnp ] = & GetMutableFextNp();

  xstate_vectors_[ state_vec_freactnp ] = & GetMutableFreactNp();

  xstate_vectors_[ state_vec_finertialnp ] = & GetMutableFinertialNp();
  xstate_vectors_[ state_vec_fviscon ] = & GetMutableFviscoN();
  xstate_vectors_[ state_vec_fvisconp ] = & GetMutableFviscoNp();

  xstate_vectors_[ state_vec_fstructold ] = & GetMutableFstructureOld();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::XStructureState::RegisterStateMatricesInMap()
{
  xstate_matrices_[ state_mat_jac ] = & GetMutableJacobian();
  xstate_matrices_[ state_mat_stiff ] = & StiffPtr();
  xstate_matrices_[ state_mat_mass ] = & GetMutableMassMatrix();
  xstate_matrices_[ state_mat_damp ] = & GetMutableDampMatrix();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::XStructureState::SetNewState( const XFEM::XFieldState & xstate )
{
  try
  {
    const XSTR::XStructureState & xstr_state =
        dynamic_cast<const XSTR::XStructureState & >( xstate );
    SetNewState( xstr_state );
  }
  catch (const std::bad_cast& e)
  {
    dserror("Dynamic cast to \"const XSTR::XStructureState\" failed!\\"
            "( throw = \" %s \" )", e.what() );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::XStructureState::SetNewState( const XSTR::XStructureState & xstate )
{
  XFEM::XFieldState::SetNewState( xstate );

  // copy vector pointers
  {
    XStateVecMap::const_iterator other_cit = xstate.XStateVectorBegin();
    for ( XStateVecMap::const_iterator cit = this->xstate_vectors_.begin();
          cit != this->xstate_vectors_.end(); ++cit )
    {
      Teuchos::RCP<Epetra_Vector> & xstate_ptr = *( cit->second );
      xstate_ptr = *( other_cit->second );

      ++other_cit;
    }
  }

  // copy matrix pointers
  {
    XStateMatMap::const_iterator other_cit = xstate.XStateMatrixBegin();
    for ( XStateMatMap::const_iterator cit = this->xstate_matrices_.begin();
          cit != this->xstate_matrices_.end(); ++cit )
    {
      Teuchos::RCP<LINALG::SparseOperator> & xstate_ptr = *( cit->second );
      xstate_ptr = *( other_cit->second );

      ++other_cit;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::XStructureState::ResetNonStandardDofs(
    const DRT::DiscretizationInterface & full_discret )
{
  for ( XStateVecMap::const_iterator cit = xstate_vectors_.begin();
      cit != xstate_vectors_.end(); ++cit )
  {
    Epetra_Vector & xstate_vec = **( cit->second );
    NaturalExtensionOfNonStandardDofValues( cit->first,
        full_discret, xstate_vec );
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::XStructureState::TransferToNewState(
    const DRT::DiscretizationInterface& new_discret,
    XFEM::XFieldState & new_xstate ) const
{
  XSTR::XStructureState * new_xstr_state =
      dynamic_cast<XSTR::XStructureState *>( &new_xstate );

  if ( not new_xstr_state )
    dserror( "NULL pointer!" );

  const Epetra_BlockMap & old_dof_row_map = (*XStateVectorBegin()->second)->Map();

  Teuchos::RCP<Epetra_Map> new_std_dof_row_map = BuildNewStandardDofRowMap(
      new_discret );

  Teuchos::RCP<Epetra_Export> old2new_dof_row_exporter =
      Teuchos::rcp( new Epetra_Export( old_dof_row_map, *new_std_dof_row_map ) );

  // state vector ( only standard dofs )
  Teuchos::RCP<Epetra_Vector> xstate_std =
      Teuchos::rcp( new Epetra_Vector( *new_std_dof_row_map ) );

  XStateVecMap::iterator it_new_state = new_xstr_state->XStateVectorBegin();
  for ( XStateVecMap::const_iterator cit = xstate_vectors_.begin();
        cit != xstate_vectors_.end(); ++cit )
  {
    const Epetra_Vector & xstate_vec = **( cit->second );
    Epetra_Vector & new_xstate_vec = **( it_new_state->second );

    if ( cit->first != it_new_state->first )
    {
      dserror( "Enumerator mismatch! ( \"%s\"[OLD] != \"%s\"[NEW] )",
          StateVectorName2String( cit->first ).c_str(),
          StateVectorName2String( it_new_state->first ).c_str() );
    }

    xstate_std->PutScalar( 0.0 );
    if ( xstate_std->Export( xstate_vec, *old2new_dof_row_exporter, Insert ) )
      dserror( "Export failed!" );

    new_xstate_vec.PutScalar( 0.0 );
    STR::AssembleVector( 0.0, new_xstate_vec, 1.0, *xstate_std );

    NaturalExtensionOfNonStandardDofValues( cit->first, new_discret, new_xstate_vec );

#if 0
    std::cout << StateVectorName2String( cit->first ) << " [OLD]" << std::endl;
    std::cout << xstate_vec << std::endl;

    std::cout << StateVectorName2String( it_new_state->first ) << " [NEW]" << std::endl;
    std::cout << new_xstate_vec << std::endl;
#endif

    // increase counter
    ++it_new_state;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::XStructureState::NaturalExtensionOfNonStandardDofValues(
    enum StateVectorName state_name,
    const DRT::DiscretizationInterface & new_discret,
    Epetra_Vector & new_xstate_vec ) const
{
  switch ( state_name )
  {
    case state_vec_disnp:
    case state_vec_velnp:
    case state_vec_accnp:
    {
      NaturalExtensionOfNonStandardDofValues( new_discret, new_xstate_vec );
      break;
    }
    default:
    {
      // do nothing for the forces for now
      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::XStructureState::NaturalExtensionOfNonStandardDofValues(
    const DRT::DiscretizationInterface & new_discret,
    Epetra_Vector & new_xstate_vec ) const
{
  const unsigned num_std_dofs = new_discret.NumStandardDof(
      0, new_discret.lRowNode( 0 ) );

  const unsigned my_num_nodes = new_discret.NodeRowMap()->NumMyElements();

  std::vector<int> dofs_std;
  std::vector<int> dofs_non_std;

  for ( unsigned lid = 0; lid < my_num_nodes; ++lid )
  {
    DRT::Node* node = new_discret.lRowNode( lid );

    const unsigned num_dofs = new_discret.NumDof( 0, node );

    // skip nodes with standard dofs only
    if ( num_dofs == num_std_dofs )
      continue;

    if ( num_dofs % num_std_dofs != 0 )
      dserror( "The dof number is no integer multiple of standard dof number!" );

    unsigned num_nodal_dofsets = num_dofs / num_std_dofs;

    dofs_std.clear();
    new_discret.Dof(dofs_std, node, 0, 1 );

    double * xstate_values = new_xstate_vec.Values();

    for ( unsigned nds = 1; nds < num_nodal_dofsets; ++nds )
    {
      dofs_non_std.clear();
      new_discret.Dof( dofs_non_std, node, 0, nds );

      for ( unsigned d = 0; d < dofs_non_std.size(); ++d )
      {
        int std_dof_gid     = dofs_std[ d ];
        int non_std_dof_gid = dofs_non_std[ d ];

        int std_dof_lid = new_xstate_vec.Map().LID( std_dof_gid );
        int non_std_dof_lid = new_xstate_vec.Map().LID( non_std_dof_gid );

        xstate_values[ non_std_dof_lid ] = xstate_values[ std_dof_lid ];
      }
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> XSTR::XStructureState::BuildNewStandardDofRowMap(
    const DRT::DiscretizationInterface & new_discret ) const
{
  const unsigned num_std_dofs = new_discret.NumStandardDof(
      0, new_discret.lRowNode( 0 ) );

  const unsigned my_num_nodes = new_discret.NodeRowMap()->NumMyElements();

  std::vector<int> my_std_row_dofs;
  my_std_row_dofs.reserve( num_std_dofs*my_num_nodes );

  std::vector<int> dofs;

  for ( unsigned lid = 0; lid < my_num_nodes; ++lid )
  {
    DRT::Node* node = new_discret.lRowNode( lid );
    dofs.clear();
    new_discret.Dof( dofs, node, 0, 0 );

    std::copy( dofs.begin(), dofs.end(), std::back_inserter( my_std_row_dofs ) );
  }

  return Teuchos::rcp( new Epetra_Map( -1, static_cast<int>( my_std_row_dofs.size() ),
      & my_std_row_dofs[0], 0, new_discret.Comm() ) );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool XSTR::XStructureState::Destroy()
{
  // delete all state matrices
  for ( XStateMatMap::const_iterator cit = xstate_matrices_.begin();
        cit != xstate_matrices_.end(); ++cit )
  {
    Teuchos::RCP<LINALG::SparseOperator> & state_mat = *( cit->second );
    try
    {
      XFEM::DestroyMatrix( state_mat, true );
    }
    catch ( const std::runtime_error & e )
    {
      std::ostringstream msg;
      msg << "The state matrix \"" << StateMatrixName2String( cit->first )
          << "\" could not be destroyed!\n\nCaught the following runtime error:\n"
          << e.what();
      dserror( msg.str().c_str() );
    }
  }

  // destroy all state Epetra_Vectors
  for ( XStateVecMap::const_iterator cit = xstate_vectors_.begin();
        cit != xstate_vectors_.end(); ++cit )
  {
    Teuchos::RCP<Epetra_Vector> & state_vec = *( cit->second );
    try
    {
      XFEM::DestroyRCPObject( state_vec, true );
    }
    catch ( const std::runtime_error & e )
    {
      std::ostringstream msg;
      msg << "The state vector \"" << StateVectorName2String( cit->first )
          << "\" could not be destroyed!\n\nCaught the following runtime error:\n"
          << e.what();
      dserror( msg.str().c_str() );
    }
  }

  // destroy Epetra_timer
  XFEM::DestroyRCPObject( GetMutableTimer(), true );

  XFEM::DestroyRCPObject( GetMutableMultiDis(), true );
  XFEM::DestroyRCPObject( GetMutableMultiVel(), true );
  XFEM::DestroyRCPObject( GetMutableMultiAcc(), true );

  // destroy state maps
  XFEM::DestroyRCPObject( GlobalProblemMapPtr(), true );

  // destroy remaining member variables
  CutWizardPtr() = Teuchos::null;
  ConditionManagerPtr() = Teuchos::null;
  XDofSetPtr() = Teuchos::null;

  return true;
}

