/*----------------------------------------------------------------------------*/
/**
\file xstr_xstructure_structure_state.cpp

\brief XStructure/Structure state handling ( combination of XFEM and std.
       discretizations )


\maintainer Michael Hiermeier

\date Jun 28, 2016

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "xstr_xstructure_structure_state.H"

#include "../drt_xfem/xfield_state_utils.H"
#include "../drt_xfem/xfem_multi_field_mapextractor.H"
#include "../drt_xfem/xfem_dofset.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_discret_xfem.H"



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XSTR::XStructureStructureState::XStructureStructureState()
{
  // intentionally left blank
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::XStructureStructureState::Setup()
{
  CheckInit();
  XStructureState::Setup();
}
