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

#include "../drt_xfem/xfield_state_utils.H"

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

  XFEM::XFieldState::issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool XSTR::XStructureState::Destroy()
{
  // delete all state matrices
  XFEM::DestroyMatrix( GetMutableJacobian(), true );
  XFEM::DestroyMatrix( GetMutableMassMatrix(), true );
  XFEM::DestroyMatrix( GetMutableDampMatrix(), true );

  // destroy Epetra_timer
  XFEM::DestroyRCPObject( GetMutableTimer(), true );

  // destroy all state TIMINT::TimIntMStep state objects
  XFEM::DestroyRCPObject( GetMutableMultiTime(), true );
  XFEM::DestroyRCPObject( GetMutableDeltaTime(), true );

  XFEM::DestroyRCPObject( GetMutableMultiDis(), true );
  XFEM::DestroyRCPObject( GetMutableMultiVel(), true );
  XFEM::DestroyRCPObject( GetMutableMultiAcc(), true );

  // destroy all state Epetra_Vectors
  XFEM::DestroyRCPObject( GetMutableDisNp(), true );
  XFEM::DestroyRCPObject( GetMutableVelNp(), true );
  XFEM::DestroyRCPObject( GetMutableAccNp(), true );

  XFEM::DestroyRCPObject( GetMutableFintN(), true );
  XFEM::DestroyRCPObject( GetMutableFintNp(), true );

  XFEM::DestroyRCPObject( GetMutableFextN(), true );
  XFEM::DestroyRCPObject( GetMutableFextNp(), true );

  XFEM::DestroyRCPObject( GetMutableFreactNp(), true );

  XFEM::DestroyRCPObject( GetMutableFinertialNp(), true );
  XFEM::DestroyRCPObject( GetMutableFinertialN(), true );

  XFEM::DestroyRCPObject( GetMutableFviscoNp(), true );
  XFEM::DestroyRCPObject( GetMutableFviscoN(), true );

  XFEM::DestroyRCPObject( GetMutableFstructureOld(), true );

  // destroy state maps
  XFEM::DestroyRCPObject( GlobalProblemMapPtr(), true );

  // destroy remaining member variables
  CutWizardPtr() = Teuchos::null;
  ConditionManagerPtr() = Teuchos::null;
  XDofSetPtr() = Teuchos::null;

  return true;
}
